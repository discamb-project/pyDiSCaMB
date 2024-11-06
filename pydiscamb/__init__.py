from ._pydiscamb import (
    __doc__,
    get_discamb_version,
    calculate_structure_factors_IAM,
    calculate_structure_factors_TAAM,
    DiscambWrapper,
    DiscambWrapperTests,
    FCalcMethod,
    get_table,
)

__all__ = [
    "__doc__",
    "get_discamb_version",
    calculate_structure_factors_IAM,
    calculate_structure_factors_TAAM,
    "DiscambWrapper",
    "FCalcMethod",
    "get_table",
]

def _get_TAAM_databank_directory() -> str:
    from pathlib import Path
    return str(Path(__file__).parent / "data")
