import pytest
from mmtbx.f_model import manager
from iotbx.phil import parse
from conftest import SKIP_IF_CCTBX_PR_NOT_MERGED

from pydiscamb.cctbx_interface import gradients_taam, from_scatterers_taam, pydiscamb_master_params
from pydiscamb.cctbx_interface.from_scatterers import CctbxFromScatterersResult
from pydiscamb.cctbx_interface.gradients import CctbxGradientsResult
from pydiscamb.discamb_wrapper import FCalcMethod

@pytest.fixture
def d_tyrosine(tyrosine):
    shaken = tyrosine.deep_copy_scatterers()
    shaken.shake_sites_in_place(rms_difference=0.1)
    shaken.shake_adp()
    shaken.shake_occupancies()
    f_obs = tyrosine.structure_factors(d_min=1.7).f_calc().as_amplitude_array()
    model = manager(f_obs=f_obs, xray_structure=shaken)
    target = model.target_functor()(compute_gradients=True)
    d_target_d_f_calc = target.d_target_d_f_calc_work()
    return f_obs, d_target_d_f_calc, tyrosine

@SKIP_IF_CCTBX_PR_NOT_MERGED
class TestCctbxObjects:

    

    def test_grads_init(self, d_tyrosine):
        f_obs, d_target_d_f_calc, tyrosine = d_tyrosine
        res = gradients_taam(
            tyrosine,
            None,
            f_obs.set(),
            d_target_d_f_calc,
            0,
            None,
            None,
        )
        assert hasattr(res, "d_target_d_site_frac")

    def test_grads_extra_params(self, d_tyrosine):
        f_obs, d_target_d_f_calc, tyrosine = d_tyrosine

        # Mock scope
        scope = parse(
            "discamb { include scope pydiscamb.cctbx_interface.pydiscamb_master_params }",
            process_includes=True,
        ).extract()
        scope.discamb.taam.unit_cell_charge = 100

        # Enable grads for a scatterer
        tyrosine.scatterers()[0].flags.set_grad_site(True)

        res1 = gradients_taam(
            tyrosine,
            None,
            f_obs.set(),
            d_target_d_f_calc,
            0,
            None,
            None,
        )
        res2 = gradients_taam(
            tyrosine,
            None,
            f_obs.set(),
            d_target_d_f_calc,
            0,
            None,
            scope,
        )

        # Checking all, but only the first will be different
        assert not pytest.approx(res1.d_target_d_site_frac()) == list(
            res2.d_target_d_site_frac()
        )

    def test_fcalc_init(self, tyrosine):
        f_obs = tyrosine.structure_factors(d_min=1.7).f_calc()
        res = from_scatterers_taam(
            tyrosine,
            f_obs.set(),
        )
        assert hasattr(res, "f_calc")

    def test_fcalc_extra_params(self, tyrosine):

        # Mock scope
        scope = parse(
            "discamb { include scope pydiscamb.cctbx_interface.pydiscamb_master_params }",
            process_includes=True,
        ).extract()
        scope.discamb.taam.unit_cell_charge = 100

        f_obs = tyrosine.structure_factors(d_min=1.7).f_calc()

        res1 = from_scatterers_taam(
            tyrosine,
            f_obs.set(),
        )
        res2 = from_scatterers_taam(
            tyrosine,
            f_obs.set(),
            extra_params=scope,
        )

        # Checking all, but only the first will be different
        assert not pytest.approx(res1.f_calc()) == list(res2.f_calc())


class TestResultsObjects:

    def test_grads_init(self, d_tyrosine):
        _, d_target_d_f_calc, tyrosine = d_tyrosine
        res = CctbxGradientsResult(tyrosine, d_target_d_f_calc.set(), d_target_d_f_calc, 0, FCalcMethod.IAM)
        assert hasattr(res, "packed")
    
    def test_fcalc_init(self, d_tyrosine):
        _, d_target_d_f_calc, tyrosine = d_tyrosine
        res = CctbxFromScatterersResult(tyrosine, d_target_d_f_calc.set(), FCalcMethod.IAM)
        assert hasattr(res, "f_calc")

    @pytest.mark.parametrize("site", [True, False])
    @pytest.mark.parametrize("adp", [False, "iso", "aniso"])
    @pytest.mark.parametrize("oc", [True, False])
    def test_grads_flags(self, d_tyrosine, site, adp, oc):
        _, d_target_d_f_calc, tyrosine = d_tyrosine
        xrs = tyrosine.deep_copy_scatterers()
        if site:
            xrs.scatterers()[0].flags.set_grad_site(True)
        if adp == "iso":
            xrs.scatterers()[0].flags.set_grad_u_iso(True)
        if adp == "aniso":
            xrs.scatterers()[0].convert_to_anisotropic(xrs.unit_cell())
            xrs.scatterers()[0].flags.set_grad_u_aniso(True)
        if oc:
            xrs.scatterers()[0].flags.set_grad_occupancy(True)

        res = CctbxGradientsResult(xrs, d_target_d_f_calc.set(), d_target_d_f_calc, 0, FCalcMethod.TAAM, algorithm="macromol")


        assert all(all(i == 0 for i in s) for s in res.d_target_d_site_frac()) ^ site
        assert all(s == 0 for s in res.d_target_d_occupancy()) ^ oc
        if adp == "iso":
            assert any(s != 0 for s in res.d_target_d_u_iso())
            assert all(all(i == 0 for i in s) for s in res.d_target_d_u_star())
        elif adp == "aniso":
            assert all(s == 0 for s in res.d_target_d_u_iso())
            assert any(any(i == 0 for i in s) for s in res.d_target_d_u_star())
        else:
            assert all(s == 0 for s in res.d_target_d_u_iso())
            assert all(all(i == 0 for i in s) for s in res.d_target_d_u_star())
    
    def test_grad_error_with_iso_aniso_flags(self, d_tyrosine):
        _, d_target_d_f_calc, tyrosine = d_tyrosine
        xrs = tyrosine.deep_copy_scatterers()
        xrs.scatterers()[0].convert_to_isotropic(xrs.unit_cell())
        xrs.scatterers()[0].flags.set_grad_u_aniso(True)

        with pytest.raises(ValueError, match="Attempted to compute aniso gradient for iso scatterer"):
            CctbxGradientsResult(xrs, d_target_d_f_calc.set(), d_target_d_f_calc, 0, FCalcMethod.TAAM, algorithm="macromol")
