import pytest
from mmtbx.f_model import manager
from iotbx.phil import parse

from pydiscamb.cctbx_interface import *

class TestCctbxGradients:

    @pytest.fixture
    def d_tyrosine(self, tyrosine):
        shaken = tyrosine.deep_copy_scatterers()
        shaken.shake_sites_in_place(rms_difference=0.1)
        shaken.shake_adp()
        shaken.shake_occupancies()
        f_obs = tyrosine.structure_factors(d_min=1.7).f_calc().as_amplitude_array()
        model = manager(f_obs=f_obs, xray_structure=shaken)
        target = model.target_functor()(compute_gradients=True)
        d_target_d_f_calc = target.d_target_d_f_calc_work()
        return f_obs, d_target_d_f_calc, tyrosine
    
    def test_init(self, d_tyrosine):
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

    def test_extra_params(self, d_tyrosine):
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
        assert not pytest.approx(res1.d_target_d_site_frac()) == list(res2.d_target_d_site_frac())


class TestCctbxFcalc:
    def test_init(self, tyrosine):
        f_obs = tyrosine.structure_factors(d_min=1.7).f_calc()
        res = from_scatterers_taam(
            tyrosine,
            f_obs.set(),
        )
        assert hasattr(res, "f_calc")

    def test_extra_params(self, tyrosine):

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
