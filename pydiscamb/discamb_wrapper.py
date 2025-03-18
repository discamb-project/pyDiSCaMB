from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union, overload

from cctbx import miller
from cctbx.array_family import flex
from cctbx.xray.structure import structure
from pydiscamb._cpp_module import (
    FCalcDerivatives,
    PythonInterface,
    TargetDerivatives,
    table_alias,
)
from pydiscamb.taam_parameters import get_default_databank


class FCalcMethod(Enum):
    IAM = 0
    TAAM = 1


class DiscambWrapper(PythonInterface):
    def __init__(self, xrs: structure, method: FCalcMethod = None, **kwargs):
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
        See :code:`DiscambWrapper._prepare_calculator_params` for details.

        Parameters
        ----------
        xrs : structure
            Object defining the symmetry, atom positions, ADPs, occupancy ect.
        method : FCalcMethod, optional
            Which structure factor model to use, by default FCalcMethod.IAM
        """
        calculator_params = self._prepare_calculator_params(xrs, method, kwargs)
        super().__init__(xrs, calculator_params)

        self._scatterer_flags = xrs.scatterer_flags()
        self._atomstr = _concat_scatterer_labels(xrs)
        self._unit_cell = xrs.unit_cell()

    @staticmethod
    def _prepare_calculator_params(
        xrs: structure, method: FCalcMethod, kwargs: Dict[str, str]
    ) -> Dict[str, str]:
        if method is None:
            method = FCalcMethod.IAM
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
        out.update(kwargs)
        # Discamb likes spaces instead of underscores
        out = {key.replace("_", " "): val for key, val in out.items()}
        return out

    def update_structure(self, xrs: structure):
        if (
            self._atomstr != _concat_scatterer_labels(xrs)
            or self._unit_cell != xrs.unit_cell()
        ):
            raise ValueError(
                "Incompatible structures. "
                "Must have same scatterers in the same order, and same unit cell"
            )
        self._scatterer_flags = xrs.scatterer_flags()
        super().update_structure(xrs)

    @classmethod
    def from_file(
        cls, filepath: Union[str, Path], method: FCalcMethod = FCalcMethod.IAM, **kwargs
    ) -> "DiscambWrapper":
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
        cls, pdb_code: str, method: FCalcMethod = FCalcMethod.IAM, **kwargs
    ) -> "DiscambWrapper":
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
    def _from_pdb_str(cls, pdb_str: str, method, **kwargs) -> "DiscambWrapper":
        import iotbx.pdb
        import mmtbx.model
        from libtbx.utils import null_out

        pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
        model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
        return cls(model.get_xray_structure(), method, **kwargs)

    ## Annotations
    @overload
    def f_calc(self, miller_set: None) -> flex.complex_double:
        ...

    @overload
    def f_calc(self, d_min: float) -> flex.complex_double:
        ...

    @overload
    def f_calc(self, miller_set: miller.set) -> miller.array:
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
    def d_target_d_params(self, d_target_d_f_calc: list) -> List[TargetDerivatives]:
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

        # Pack selected gradients
        size = 0
        for s in self._scatterer_flags:
            size += 3 * s.grad_site()
            size += 1 * s.grad_u_iso()
            size += 6 * s.grad_u_aniso()
            size += 1 * s.grad_occupancy()
        out = flex.double(size, 0)
        o_ind = 0
        for s_ind, s in enumerate(self._scatterer_flags):
            if s.grad_site():
                out[o_ind] = res[s_ind].site_derivatives[0]
                o_ind += 1
                out[o_ind] = res[s_ind].site_derivatives[1]
                o_ind += 1
                out[o_ind] = res[s_ind].site_derivatives[2]
                o_ind += 1
            if s.grad_u_iso():
                out[o_ind] = res[s_ind].adp_derivatives[0]
                o_ind += 1
            if s.grad_u_aniso():
                out[o_ind] = res[s_ind].adp_derivatives[0]
                o_ind += 1
                out[o_ind] = res[s_ind].adp_derivatives[1]
                o_ind += 1
                out[o_ind] = res[s_ind].adp_derivatives[2]
                o_ind += 1
                out[o_ind] = res[s_ind].adp_derivatives[3]
                o_ind += 1
                out[o_ind] = res[s_ind].adp_derivatives[4]
                o_ind += 1
                out[o_ind] = res[s_ind].adp_derivatives[5]
                o_ind += 1
            if s.grad_occupancy():
                out[o_ind] = res[s_ind].occupancy_derivatives
                o_ind += 1
        assert o_ind == size
        return out


class DiscambWrapperCached(DiscambWrapper):

    __cache = {}

    def __init__(self, xrs: structure, method: FCalcMethod = None, **kwargs):
        if self.__check_cache(xrs, method, kwargs) is not None:
            self.update_structure(xrs)
            return

        super().__init__(xrs, method, **kwargs)

        self.__cache[self.__get_cache_key(xrs, method, kwargs)] = self

    def __new__(cls, xrs: structure, method: FCalcMethod = None, **kwargs):
        cache = cls.__check_cache(xrs, method, kwargs)
        if cache is None:
            return super().__new__(cls)
        return cache

    @classmethod
    def __check_cache(cls, xrs: structure, method: FCalcMethod, kwargs: Dict[str, str]):
        key = cls.__get_cache_key(xrs, method, kwargs)
        return cls.__cache.get(key)

    @classmethod
    def __get_cache_key(
        cls, xrs: structure, method: FCalcMethod, kwargs: Dict[str, str]
    ):
        atomstr = _concat_scatterer_labels(xrs)
        unitcell = xrs.unit_cell()
        params = cls._prepare_calculator_params(xrs, method, kwargs)
        key = (atomstr, unitcell, *sorted(params.items()))
        return key


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


def _concat_scatterer_labels(xrs: structure) -> str:
    return "".join(s.label for s in xrs.scatterers())
