from pathlib import Path
from typing import List

_DEFAULT_TAAM_DATABANK: str = None
_TAAM_ROOT: Path = None


def get_TAAM_root() -> Path:
    global _TAAM_ROOT
    if _TAAM_ROOT is None:
        # Check if editable
        project_root = Path(__file__).parent.parent
        readme = project_root / "README.md"
        if readme.exists() and readme.read_text()[:11] == "# pydiscamb":
            _TAAM_ROOT = project_root / "data"
        else:
            _TAAM_ROOT = Path(__file__).parent / "data"
        assert _TAAM_ROOT.exists()
    return _TAAM_ROOT


def get_TAAM_databanks() -> List[str]:
    """
    Get a list of all available databanks.
    """
    # Upon installation, all *databank.txt in the data-folder
    # are copied into the installation directory of the module,
    # without preserving other folder structure.
    # Assume no filenames are duplicate.
    files = get_TAAM_root().glob("**/*databank.txt")
    return [str(file) for file in files]


def is_MATTS_installed() -> bool:
    return any("MATTS" in path for path in get_TAAM_databanks())


def get_default_databank() -> str:
    global _DEFAULT_TAAM_DATABANK
    if _DEFAULT_TAAM_DATABANK is None:
        banks = get_TAAM_databanks()
        search = "MATTS" if is_MATTS_installed() else "default"
        for bank in banks:
            if search in bank:
                _DEFAULT_TAAM_DATABANK = bank
                break
        assert Path(_DEFAULT_TAAM_DATABANK).exists()
    return _DEFAULT_TAAM_DATABANK
