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

@pytest.mark.parametrize(
    "hkl",
    [
        (0, 0, 0),
        (1, 2, 3),
        (2, -2, 0),
    ]
)
def test_single_hkl(random_structure, hkl):
    w = DiscambWrapper(random_structure)
    grad_tuple = w.d_f_calc_hkl_d_params(hkl)
    h, k, l = hkl
    grad_int = w.d_f_calc_hkl_d_params(h, k, l)
    assert grad_tuple is not grad_int
    assert grad_int.hkl == grad_tuple.hkl
    assert grad_int.structure_factor == grad_tuple.structure_factor
    assert grad_int.adp_derivatives == grad_tuple.adp_derivatives
    assert grad_int.site_derivatives == grad_tuple.site_derivatives
    assert grad_int.occupancy_derivatives == grad_tuple.occupancy_derivatives
