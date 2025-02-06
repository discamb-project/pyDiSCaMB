from pathlib import Path
from typing import List

ROOT = Path(__file__).parent / "data"


def get_TAAM_databanks() -> List[str]:
    """
    Get a list of all available databanks.
    """
    # Upon installation, all *databank.txt in the data-folder
    # are copied into the installation directory of the module,
    # without preserving other folder structure.
    # Assume no filenames are duplicate.
    files = ROOT.glob("*databank.txt")
    return [str(file) for file in files]


def get_TAAM_root() -> str:
    return str(ROOT)


def is_MATTS_installed() -> bool:
    return any("MATTS" in path for path in get_TAAM_databanks())


def get_default_databank() -> str:
    banks = get_TAAM_databanks()
    search = "MATTS" if is_MATTS_installed() else "default"
    for bank in banks:
        if search in bank:
            return bank
    # Failsafe
    return bank
