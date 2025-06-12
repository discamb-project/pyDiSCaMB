from cctbx.xray.structure_factors.from_scatterers_direct import from_scatterers_direct
from cctbx.xray.structure_factors.manager import managed_calculation_base

from pydiscamb.cctbx_interface.phil_scope import scope_to_taam_dict
from pydiscamb.discamb_wrapper import DiscambWrapperCached, FCalcMethod


class from_scatterers_taam(from_scatterers_direct):

    def __init__(
        self,
        xray_structure,
        miller_set,
        manager=None,
        algorithm="taam",
        extra_params=None,
    ):
        if extra_params is None:
            extra_params_dict = {}
        else:
            extra_params_dict = scope_to_taam_dict(extra_params)

        # TODO add timings
        managed_calculation_base.__init__(
            self, manager, xray_structure, miller_set, algorithm="taam"
        )
        self._results = CctbxFromScatterersResult(
            xray_structure, miller_set, FCalcMethod.TAAM, **extra_params_dict
        )


class CctbxFromScatterersResult:
    def __init__(self, xrs, miller_set, method, **kwargs):
        w = DiscambWrapperCached(xrs, method, **kwargs)
        w.set_indices(miller_set.indices())
        self._fcalc = w.f_calc()

    def f_calc(self):
        return self._fcalc
