import pytest

from pydiscamb import DiscambWrapper

def test_f_calc(random_structure):
    w = DiscambWrapper(random_structure)
    indices = [
        (0, 1, 2),
        (2, 3, 4),
        (10, 11, 12)
    ]
    w.set_indices(indices)

    fcalc = w.f_calc()
    grad = w.d_f_calc_d_params()

    assert pytest.approx(fcalc) == [g.structure_factor for g in grad]

def test_indices(random_structure):
    w = DiscambWrapper(random_structure)
    indices = [
        (0, 1, 2),
        (2, 3, 4),
        (10, 11, 12)
    ]
    w.set_indices(indices)
    grad = w.d_f_calc_d_params()


    assert indices == [tuple(g.hkl) for g in grad]

def test_output_sizes(random_structure):
    w = DiscambWrapper(random_structure)
    indices = [
        (0, 1, 2),
        (2, 3, 4),
        (10, 11, 12)
    ]
    w.set_indices(indices)
    grad = w.d_f_calc_d_params()

    n_scatterers = random_structure.scatterers().size()
    n_adps = 1 if random_structure.use_u_iso() else 6

    assert len(grad) == len(indices)
    assert all(len(g.site_derivatives) == n_scatterers for g in grad)
    assert all(len(g.adp_derivatives) == n_scatterers * n_adps for g in grad)
    assert all(len(g.occupancy_derivatives) == n_scatterers for g in grad)

def test_simple_structure_site_only():
    from cctbx.development import random_structure as cctbx_random_structure
    from cctbx.sgtbx import space_group_info

    group = space_group_info(19)
    xrs = cctbx_random_structure.xray_structure(
        space_group_info=group,
        elements=["Au"],
        general_positions_only=False,
        use_u_iso=True,
        random_u_iso=False,
        random_occupancy=False,
    )
    xrs.scattering_type_registry(table="electron")
    xrs.scatterers()[0].site = (0, 0, 0)

    w = DiscambWrapper(xrs)
    w.set_d_min(2)
    f_calc = w.f_calc()

    def target(fobs, fcalc):
        return sum((fo - fc)**2 for fo, fc in zip(fobs, fcalc))
    
    def d_target_d_fcalc(fobs, fcalc):
        return [2*fc - 2*fo for fo, fc in zip(fobs, fcalc)]
    
    def d_target_d_params(fobs)
    grad = w.d_f_calc_d_params()