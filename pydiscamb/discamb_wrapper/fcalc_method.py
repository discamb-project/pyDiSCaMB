import os
import tempfile
from enum import Enum
from typing import TYPE_CHECKING, Dict, Tuple, List


from pydiscamb._cpp_module import table_alias
from pydiscamb.taam_parameters import get_default_databank

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
                "algorithm": "standard",
            }
        elif self == FCalcMethod.TAAM:
            out = {
                "model": "taam",
                "bank_path": get_default_databank(),
                "electron_scattering": table == "electron",
                "table": alias,
                "assignment_csv": _get_tmp_assignment_filename(),
                "algorithm": "macromol",
            }

        out.update(kwargs)
        # Discamb likes spaces instead of underscores
        out = {key.replace("_", " "): val for key, val in out.items()}
        return out

    def to_cache_lookup_key(
        self, xrs: "structure", kwargs: Dict[str, str]
    ) -> List[Tuple[str, str]]:
        dct = self.to_dict(xrs, kwargs)
        if self == FCalcMethod.IAM:
            # No issues with IAM
            pass
        elif self == FCalcMethod.TAAM:
            # We might have a new assignment csv being generated
            if kwargs.get("assignment_csv") is None:
                try:
                    os.remove(dct.pop("assignment csv"))
                except PermissionError:
                    # Known failure on windows
                    pass

        return sorted(dct.items())
