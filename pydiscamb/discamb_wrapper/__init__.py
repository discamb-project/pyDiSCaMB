from typing import List

from pydiscamb.discamb_wrapper.discamb_wrapper import DiscambWrapper
from pydiscamb.discamb_wrapper.cache import DiscambWrapperCached
from pydiscamb.discamb_wrapper.fcalc_method import FCalcMethod

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cctbx.xray.structure import structure


__all__ = [
    "DiscambWrapper",
    "DiscambWrapperCached",
    "FCalcMethod",
    "calculate_structure_factors_IAM",
    "calculate_structure_factors_TAAM",
]


def _calculate_structure_factors(
    xrs: "structure", d_min: float, method: FCalcMethod
) -> List[complex]:
    w = DiscambWrapper(xrs, method)
    w.set_d_min(d_min)
    return w.f_calc()


def calculate_structure_factors_IAM(xrs: "structure", d_min: float) -> List[complex]:
    return _calculate_structure_factors(xrs, d_min, FCalcMethod.IAM)


def calculate_structure_factors_TAAM(xrs: "structure", d_min: float) -> List[complex]:
    return _calculate_structure_factors(xrs, d_min, FCalcMethod.TAAM)
