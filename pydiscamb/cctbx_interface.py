from __future__ import absolute_import, division, print_function
from cctbx.xray.structure_factors.gradients_base import gradients_base
from cctbx.xray.structure_factors.misc import expensive_function_call_message
from cctbx import adptbx
from libtbx.utils import user_plus_sys_time
from cctbx.array_family import flex

from pydiscamb.discamb_wrapper import DiscambWrapper, FCalcMethod


class gradients_discamb(gradients_base):

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
            algorithm="discamb",
        )
        w = DiscambWrapper(xray_structure, FCalcMethod.TAAM)
        w.set_indices(miller_set.indices())
        self._grads = w.d_target_d_params(list(d_target_d_f_calc))

        self.d_target_d_site_cart_was_used = False
        self.d_target_d_u_cart_was_used = False

    def d_target_d_site_frac(self):
        return self.check_size(
            flex.vec3_double([g.site_derivatives for g in self._grads])
        )

    def d_target_d_u_iso(self):
        return self.check_size(
            flex.double(
                [
                    g.site_derivatives[0] if len(g.site_derivatives) == 1 else 0
                    for g in self._grads
                ]
            )
        )

    def d_target_d_occupancy(self):
        return self.check_size(
            flex.double([g.occupancy_derivatives for g in self._grads])
        )

    def d_target_d_fp(self):
        return self.check_size(flex.double([0 for g in self._grads]))

    def d_target_d_fdp(self):
        return self.check_size(flex.double([0 for g in self._grads]))

    def packed(self):
        size = 0
        for s in self.xray_structure().scatterer_flags():
            size += 3 * s.grad_site()
            size += 1 * s.grad_u_iso()
            size += 6 * s.grad_u_aniso()
            size += 1 * s.grad_occupancy()
        out = flex.double(size, 0)
        o_ind = 0
        for s_ind, s in enumerate(self.xray_structure().scatterer_flags()):
            if s.grad_site():
                for j in range(3):
                    out[o_ind] = self._grads[s_ind].site_derivatives[j]
                    o_ind += 1
            if s.grad_u_iso():
                out[o_ind] = self._grads[s_ind].adp_derivatives[0]
                o_ind += 1
            if s.grad_u_aniso():
                for j in range(6):
                    out[o_ind] = self._grads[s_ind].adp_derivatives[j]
                    o_ind += 1
            if s.grad_occupancy():
                out[o_ind] = self._grads[s_ind].occupancy_derivatives
                o_ind += 1
        assert o_ind == size
        return out

    def d_target_d_site_cart(self):
        if self.d_target_d_site_cart_was_used:
            raise RuntimeError(expensive_function_call_message)
        self.d_target_d_site_cart_was_used = True
        return (
            self.d_target_d_site_frac()
            * self.xray_structure().unit_cell().fractionalization_matrix()
        )

    def d_target_d_u_star(self):
        return self.check_size(
            flex.sym_mat3(
                [
                    g.site_derivatives if len(g.site_derivatives) == 6 else [0] * 6
                    for g in self._grads
                ]
            )
        )

    def d_target_d_u_cart(self):
        if self.d_target_d_u_cart_was_used:
            raise RuntimeError(expensive_function_call_message)
        self.d_target_d_u_cart_was_used = True
        return adptbx.grad_u_star_as_u_cart(
            self.xray_structure().unit_cell(), self.d_target_d_u_star()
        )
