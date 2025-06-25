from pathlib import Path
from typing import TypeVar, Union

from pydiscamb.discamb_wrapper.fcalc_method import FCalcMethod

Self = TypeVar("Self", bound="FactoryMethodsMixin")


class FactoryMethodsMixin:

    @classmethod
    def from_file(
        cls: type[Self],
        filepath: Union[str, Path],
        method: FCalcMethod = FCalcMethod.IAM,
        **kwargs,
    ) -> Self:
        """Read a structure from a file

        Parameters
        ----------
        pdb_file : Union[str, Path]
            Path to structure to read
        method : FCalcMethod, optional
            How to calculate structure factors, by default FCalcMethod.IAM
        **kwargs :
            key-value parameters sent to `discamb/Scattering/SfCalculator::create`

        Returns
        -------
        DiscambWrapper
            Using the structure read from file

        Raises
        ------
        FileNotFoundError
            If given file is not found
        ValueError
            If unsupported file format if given
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError
        if filepath.suffix.lower() not in [".cif", ".mmcif", ".pdb"]:
            raise ValueError(
                f"Supported files are .cif, .mmcif and .pdb. Got {filepath.suffix}"
            )

        # Try reading as "normal" cif, might be mmcif
        if filepath.suffix == ".cif":
            import iotbx.cif
            cif = iotbx.cif.reader(input_string=filepath.read_text())
            structures = cif.build_crystal_structures()
            if len(structures.keys()) == 1:
                return cls(structures.popitem()[1], method=method, **kwargs)
            if len(structures.keys()) > 1:
                raise ValueError("Multiple structures found in file")
            # len is 0, assume mmcif

        return cls._from_pdb_str(filepath.read_text(), method=method, **kwargs)

    @classmethod
    def from_pdb_code(
        cls: type[Self], pdb_code: str, method: FCalcMethod = FCalcMethod.IAM, **kwargs
    ) -> Self:
        """Download structure from rcsb.org

        Parameters
        ----------
        pdb_code : str
            PDB code on rcsb.org
        method : FCalcMethod, optional
            How to calculate structure factors, by default FCalcMethod.IAM
        **kwargs :
            key-value parameters sent to `discamb/Scattering/SfCalculator::create`

        Returns
        -------
        DiscambWrapper
            Wrapper object from the structure on

        Raises
        ------
        ValueError
            _description_
        """
        if len(pdb_code) != 4:
            raise ValueError("pdb code must be 4 characters long")
        if not pdb_code.isalnum():
            raise ValueError("pdb code must be alphanumeric")

        import requests

        response = requests.get(f"https://files.rcsb.org/view/{pdb_code}.pdb")
        if response.status_code == 404:
            raise ValueError("pdb code not found on rcsb.org")
        elif response.status_code != 200:
            raise RuntimeError("Communication error with rcsb.org")
        pdb_str = response.content.decode("utf-8")
        return cls._from_pdb_str(pdb_str, method, **kwargs)

    @classmethod
    def _from_pdb_str(cls: type[Self], pdb_str: str, method, **kwargs) -> Self:
        import iotbx.pdb
        import mmtbx.model
        from libtbx.utils import null_out
        pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
        model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
        return cls(model.get_xray_structure(), method, **kwargs)
