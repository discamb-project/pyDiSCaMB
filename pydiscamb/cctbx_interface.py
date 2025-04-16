from cctbx.xray.structure_factors.gradients_base import gradients_base
from cctbx.xray.structure_factors.gradients_direct import gradients_direct
from cctbx.xray.structure_factors.manager import managed_calculation_base
from cctbx.xray.structure_factors.from_scatterers_direct import from_scatterers_direct
from cctbx.array_family import flex

from pydiscamb.discamb_wrapper import DiscambWrapperCached, FCalcMethod


class gradients_taam(gradients_direct):
    # TODO consider moving this class to cctbx
    # TODO add timings
    def __init__(
        self,
        xray_structure,
        u_iso_refinable_params,
        miller_set,
        d_target_d_f_calc,
        n_parameters,
        manager=None,
        cos_sin_table=False,
    ):
        gradients_base.__init__(
            self,
            manager,
            xray_structure,
            miller_set,
            algorithm="taam",
        )
        self._results = CctbxGradientsResult(
            self.xray_structure(),
            miller_set,
            d_target_d_f_calc,
            n_parameters,
            FCalcMethod.TAAM,
        )
        self.d_target_d_site_cart_was_used = False
        self.d_target_d_u_cart_was_used = False


class from_scatterers_taam(from_scatterers_direct):

    def __init__(
        self,
        xray_structure,
        miller_set,
        manager=None,
        cos_sin_table=False,
        algorithm="taam",
    ):
        # TODO add timings
        managed_calculation_base.__init__(
            self, manager, xray_structure, miller_set, algorithm="taam"
        )
        self._results = CctbxFromScatterersResult(
            xray_structure, miller_set, FCalcMethod.TAAM
        )


class CctbxGradientsResult:
    def __init__(
        self, xrs, miller_set, d_target_d_f_calc, n_parameters, method, **kwargs
    ):
        w = DiscambWrapperCached(xrs, method, **kwargs)
        w.set_indices(miller_set.indices())
        grads = w.d_target_d_params(list(d_target_d_f_calc))

        self._d_target_d_site_frac = flex.vec3_double(
            xrs.scatterers().size(), (0, 0, 0)
        )
        self._d_target_d_u_iso = flex.double(xrs.scatterers().size(), 0)
        self._d_target_d_u_star = flex.sym_mat3_double(
            xrs.scatterers().size(), (0, 0, 0, 0, 0, 0)
        )
        self._d_target_d_occupancy = flex.double(xrs.scatterers().size(), 0)
        self._d_target_d_fp = flex.double(xrs.scatterers().size(), 0)
        self._d_target_d_fdp = flex.double(xrs.scatterers().size(), 0)

        for i, s in enumerate(xrs.scatterer_flags()):
            if s.grad_site():
                self._d_target_d_site_frac[i] = grads[i].site_derivatives
            if s.grad_u_iso():
                self._d_target_d_u_iso[i] = grads[i].adp_derivatives[0]
            if s.grad_u_aniso():
                for j in range(6):
                    self._d_target_d_u_star[i] = grads[i].adp_derivatives[j]
            if s.grad_occupancy():
                self._d_target_d_occupancy[i] = grads[i].occupancy_derivatives

        self._packed = flex.double(n_parameters, 0)
        if n_parameters <= 0:
            return

        packed_ind = 0
        for i, s in enumerate(xrs.scatterer_flags()):
            if s.grad_site():
                for j in range(3):
                    self._packed[packed_ind] = grads[i].site_derivatives[j]
                    packed_ind += 1
            if s.grad_u_iso():
                self._packed[packed_ind] = grads[i].adp_derivatives[0]
                packed_ind += 1
            if s.grad_u_aniso():
                for j in range(6):
                    self._packed[packed_ind] = grads[i].adp_derivatives[j]
                    packed_ind += 1
            if s.grad_occupancy():
                self._packed[packed_ind] = grads[i].occupancy_derivatives
                packed_ind += 1
        assert packed_ind == n_parameters

    def d_target_d_site_frac(self):
        return self._d_target_d_site_frac

    def d_target_d_u_iso(self):
        return self._d_target_d_u_iso

    def d_target_d_u_star(self):
        return self._d_target_d_u_star

    def d_target_d_occupancy(self):
        return self._d_target_d_occupancy

    def d_target_d_fp(self):
        return self._d_target_d_fp

    def d_target_d_fdp(self):
        return self._d_target_d_fdp

    def packed(self):
        return self._packed


class CctbxFromScatterersResult:
    def __init__(self, xrs, miller_set, method, **kwargs):
        w = DiscambWrapperCached(xrs, method, **kwargs)
        w.set_indices(miller_set.indices())
        self._fcalc = w.f_calc()

    def f_calc(self):
        return self._fcalc
