from ._pydiscamb import (
    __doc__,
    get_discamb_version,
    calculate_structure_factors_IAM,
    calculate_structure_factors_TAAM,
    DiscambWrapper,
    DiscambWrapperTests,
    FCalcMethod,
)

__all__ = [
    "__doc__",
    "get_discamb_version",
    calculate_structure_factors_IAM,
    calculate_structure_factors_TAAM,
    "DiscambWrapper",
    "FCalcMethod",
]

def _get_TAAM_databank_directory() -> str:
    from pathlib import Path
    return str(Path(__file__).parent / "data")
