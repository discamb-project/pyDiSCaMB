import pytest

from pydiscamb import DiscambWrapper, FCalcMethod

@pytest.mark.parametrize("method", [FCalcMethod.IAM, FCalcMethod.TAAM])
def test_f_calc(tyrosine, method):
    w = DiscambWrapper(tyrosine, method=method)
    indices = [
        (0, 1, 2),
        (2, 3, 4),
        (10, 11, 12)
    ]
    w.set_indices(indices)

    fcalc = w.f_calc()
    grad = w.d_f_calc_d_params()

    assert pytest.approx(fcalc) == [g.structure_factor for g in grad]

@pytest.mark.parametrize("method", [FCalcMethod.IAM, FCalcMethod.TAAM])
def test_indices(tyrosine, method):
    w = DiscambWrapper(tyrosine, method=method)
    indices = [
        (0, 1, 2),
        (2, 3, 4),
        (10, 11, 12)
    ]
    w.set_indices(indices)
    grad = w.d_f_calc_d_params()


    assert indices == [tuple(g.hkl) for g in grad]

@pytest.mark.parametrize("method", [FCalcMethod.IAM, FCalcMethod.TAAM])
def test_output_sizes(tyrosine, method):
    w = DiscambWrapper(tyrosine, method=method)
    indices = [
        (0, 1, 2),
        (2, 3, 4),
        (10, 11, 12)
    ]
    w.set_indices(indices)
    grad = w.d_f_calc_d_params()

    n_scatterers = tyrosine.scatterers().size()
    n_adps = 1 if tyrosine.use_u_iso() else 6

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
        (10, 20, 30),
    ]
)
@pytest.mark.parametrize("method", [FCalcMethod.IAM, FCalcMethod.TAAM])
def test_single_hkl(tyrosine, hkl, method):
    w = DiscambWrapper(tyrosine, method=method)
    grad_tuple = w.d_f_calc_hkl_d_params(hkl)
    h, k, l = hkl
    grad_int = w.d_f_calc_hkl_d_params(h, k, l)
    assert grad_tuple is not grad_int
    assert grad_int.hkl == grad_tuple.hkl
    assert grad_int.structure_factor == grad_tuple.structure_factor
    assert grad_int.adp_derivatives == grad_tuple.adp_derivatives
    assert grad_int.site_derivatives == grad_tuple.site_derivatives
    assert grad_int.occupancy_derivatives == grad_tuple.occupancy_derivatives

@pytest.mark.parametrize(
    "hkl",
    [
        (0, 0, 1),
        (1, 2, 3),
        (2, -2, 0),
        (10, 20, 30),
    ]
)
def test_iam_vs_taam(tyrosine, hkl):
    iam = DiscambWrapper(tyrosine, method=FCalcMethod.IAM).d_f_calc_hkl_d_params(hkl)
    taam = DiscambWrapper(tyrosine, method=FCalcMethod.TAAM).d_f_calc_hkl_d_params(hkl)

    assert iam.hkl == taam.hkl
    assert pytest.approx(iam.structure_factor) != taam.structure_factor
    assert pytest.approx([i[0] for i in iam.adp_derivatives]) != [i[0] for i in taam.adp_derivatives]
    
    h, k, l = hkl
    if h:
        assert pytest.approx([i[0] for i in iam.site_derivatives]) != [i[0] for i in taam.site_derivatives]
    if k:
        assert pytest.approx([i[1] for i in iam.site_derivatives]) != [i[1] for i in taam.site_derivatives]
    if l:
        assert pytest.approx([i[2] for i in iam.site_derivatives]) != [i[2] for i in taam.site_derivatives]
    assert pytest.approx(iam.occupancy_derivatives) != taam.occupancy_derivatives
