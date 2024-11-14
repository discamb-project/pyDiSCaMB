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
