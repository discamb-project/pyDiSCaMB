import tempfile
from enum import Enum

from pydiscamb._cpp_module import (
    table_alias,
)
from pydiscamb.taam_parameters import get_default_databank

from typing import TYPE_CHECKING, Dict

if TYPE_CHECKING:
    from cctbx.xray.structure import structure


def _get_tmp_assignment_filename() -> str:
    _, out = tempfile.mkstemp(".csv", "pyDiSCaMB-atom-assignment")
    return out


class FCalcMethod(Enum):
    IAM = 0
    TAAM = 1

    def to_dict(self, xrs: "structure", kwargs: Dict[str, str]) -> Dict[str, str]:
        table = xrs.get_scattering_table()
        alias = table_alias("xray" if table is None else table)
        if self == FCalcMethod.IAM:
            out = {
                "model": "iam",
                "table": alias,
                "electron_scattering": False,
            }
        elif self == FCalcMethod.TAAM:
            out = {
                "model": "taam",
                "bank_path": get_default_databank(),
                "electron_scattering": table == "electron",
                "table": alias,
                "assignment_csv": _get_tmp_assignment_filename(),
            }

        out.update(kwargs)
        # Discamb likes spaces instead of underscores
        out = {key.replace("_", " "): val for key, val in out.items()}
        return out
