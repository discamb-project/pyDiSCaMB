from typing import overload, Tuple, Dict, Any, List
from enum import Enum

from cctbx.xray.structure import structure
from cctbx.array_family import flex
from cctbx import miller

from pydiscamb._cpp_module import PythonInterface, FCalcDerivatives, table_alias
from pydiscamb.taam_parameters import get_default_databank


class FCalcMethod(Enum):
    IAM = 0
    TAAM = 1


def get_default_calculator_params(
    xrs: structure, method: FCalcMethod
) -> Dict[str, Any]:
    table = xrs.get_scattering_table()
    out = {
        "model": "iam",
        "bank_path": get_default_databank(),
        "table": table_alias("xray" if table is None else table),
        "electron_scattering": table == "electron",
    }
    if method == FCalcMethod.IAM:
        out.update({"electron_scattering": False})
    elif method == FCalcMethod.TAAM:
        out.update({"model": "taam"})
    return out


class DiscambWrapper(PythonInterface):

    def __init__(self, xrs: structure, method=FCalcMethod.IAM, **kwargs):
        """
        Initialize a wrapper object for structure factor calculations using DiSCaMB.
        Pass an xray structure and either a FCalcMethod, 
        or kwargs which are passed to  C++: :code:`discamb::SfCalculator::create`,
        found in discamb/Scattering/SfCalculator.h

        Notes
        -----
        The :code:`method` parameter is used for lookup of sensible default
        values to be sent to :code:`discamb::SfCalculator::create`.
        Any kwargs present will override the defaults.
        The defaults are e.g. inferring the scattering table from `xrs`,
        and using the default TAAM databank.
        See :code:`get_default_calculator_params` for details.

        Parameters
        ----------
        xrs : structure
            Object defining the symmetry, atom positions, ADPs, occupancy ect.
        method : FCalcMethod, optional
            Which structure factor model to use, by default FCalcMethod.IAM
        """
        calculator_params = get_default_calculator_params(xrs, method)
        calculator_params.update(kwargs)
        # Discamb likes spaces instead of underscores
        calculator_params = {
            key.replace("_", " "): val for key, val in calculator_params.items()
        }

        super().__init__(xrs, calculator_params)
        self._scatterer_flags = xrs.scatterer_flags()

    ## Annotations
    @overload
    def f_calc(self, miller_array: None) -> flex.complex_double:
        ...

    @overload
    def f_calc(self, d_min: float) -> flex.complex_double:
        ...

    @overload
    def f_calc(self, miller_array: miller.set) -> miller.array:
        ...

    @overload
    def d_f_calc_hkl_d_params(self, hkl: Tuple[int, int, int]) -> FCalcDerivatives:
        ...

    @overload
    def d_f_calc_hkl_d_params(self, h: int, k: int, l: int) -> FCalcDerivatives:
        ...

    @overload
    def d_target_d_params(self, d_target_d_f_calc: miller.array) -> flex.double:
        ...

    @overload
    def d_target_d_params(self, d_target_d_f_calc: list) -> List[FCalcDerivatives]:
        ...

    ## Implementations
    def f_calc(self, miller_set=None):
        if miller_set is None:
            return flex.complex_double(super().f_calc())
        if isinstance(miller_set, (int, float)):  # d_min
            self.set_d_min(miller_set)
            return self.f_calc()
        if not isinstance(miller_set, miller.set):
            raise ValueError("`miller_set` must be of type `cctbx.miller.set")

        self.set_indices(miller_set.indices())
        res = self.f_calc()
        return miller_set.array(data=res)

    def d_f_calc_hkl_d_params(self, h, k=None, l=None):
        if isinstance(h, tuple):
            assert len(h) == 3
            h, k, l = h
            return self.d_f_calc_hkl_d_params(h, k, l)
        return super().d_f_calc_hkl_d_params(h, k, l)

    def d_target_d_params(self, d_target_d_f_calc):
        if isinstance(d_target_d_f_calc, list):
            return super().d_target_d_params(d_target_d_f_calc)

        self.set_indices(d_target_d_f_calc.indices())
        res = self.d_target_d_params(list(d_target_d_f_calc.data()))

        # Select desired gradients
        s = self._scatterer_flags[0]
        if s.grad_site():
            return flex.vec3_double([g.site_derivatives for g in res]).as_double()
        if s.grad_u_iso():
            return flex.double([g.adp_derivatives[0] for g in res])

        raise ValueError("Gradient flags not supported")


def _calculate_structure_factors(
    xrs: structure, d_min: float, method: FCalcMethod
) -> List[complex]:
    w = DiscambWrapper(xrs, method)
    w.set_d_min(d_min)
    return w.f_calc()


def calculate_structure_factors_IAM(xrs: structure, d_min: float) -> List[complex]:
    return _calculate_structure_factors(xrs, d_min, FCalcMethod.IAM)


def calculate_structure_factors_TAAM(xrs: structure, d_min: float) -> List[complex]:
    return _calculate_structure_factors(xrs, d_min, FCalcMethod.TAAM)
