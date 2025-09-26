from time import perf_counter, sleep
from pathlib import Path
from functools import reduce, wraps
from typing import Callable, Union
import dataclasses

from cctbx.array_family import flex
from cctbx.xray import structure, miller, crystal
from mmtbx.f_model import manager as f_model
from tqdm import trange, tqdm
from cctbx.development import random_structure as cctbx_random_structure
import numpy as np
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out

from pydiscamb import DiscambWrapper, FCalcMethod
from cctbx.xray import structure_factors


def measure(func: Callable) -> Callable:
    """Decorator for member functions of `RuntimeBase` to measure the runtime"""
    funcname = func.__name__

    @wraps(func)
    def inner(self: RuntimeBase, *args, **kwargs):
        start = perf_counter()
        out = func(self, *args, **kwargs)
        end = perf_counter()
        self.runtimes.append(
            RuntimeResult(
                start,
                end,
                self.__class__.__name__,
                funcname,
                self.xrs.scatterers().size(),
                self.miller_set.size(),
            )
        )
        return out

    return inner


@dataclasses.dataclass
class RuntimeResult:
    """Store the result of a runtime measurement"""

    start: float
    end: float
    cls: str
    method: str
    n_scatterers: int
    n_reflection: int
    time: float = dataclasses.field(init=False)

    def __post_init__(self):
        self.time = self.end - self.start

    def csv(self) -> str:
        return ",".join(
            str(getattr(self, field.name)) for field in dataclasses.fields(self)
        )

    @classmethod
    def csv_header(cls) -> str:
        return ",".join(field.name for field in dataclasses.fields(cls))


CALCS: list[type["RuntimeBase"]] = []
"""List of all subclasses of `RuntimeBase`"""


class RuntimeBase:

    def __init__(
        self, xrs: structure, d_min: float = None, miller_set: miller.set = None
    ):
        """Abstract class to measure the runtime of computing structure factors and gradients for a given structure and miller set

        Parameters
        ----------
        xrs : structure
            Structure to compute for
        d_min : float, optional
            Resolution, used to generate miller indices. If None, `miller_set` must be given, by default None
        miller_set : miller.set, optional
            Reflections to compute for. If None, `d_min` must be set, by default None
        """
        assert d_min is not None or miller_set is not None
        self.xrs = xrs
        if d_min is not None:
            miller_set = self.xrs.build_miller_set(False, d_min=d_min)
        self.miller_set = miller_set
        self.runtimes: list[RuntimeResult] = []

    def __init_subclass__(cls):
        # Register subclasses
        CALCS.append(cls)

    def f_calc(self) -> None:
        raise NotImplementedError

    def grads(self, d_target_d_f_calc: flex.complex_double) -> None:
        raise NotImplementedError

    def run(
        self,
        n: int = 1,
        show_pbar: bool = False,
        pbar_position: int = 0,
        pbar_leave: bool = False,
        overall_pbar: tqdm = None,
    ) -> None:
        """Run `self.f_calc` and `self.grads` a specified number of times

        Parameters
        ----------
        n : int, optional
            Number of runs, by default 1
        show_pbar : bool
            Whether to show a progress bar, by default False
        pbar_position : int, optional
            Sent to tqdm (position=pbar_position), by default 0
        pbar_leave : bool, optional
            Sent to tqdm (leave=pbar_leave), by default False
        overall_pbar : tqdm, optional
            Initialized tqdm object, will run `overall_pbar.update()` on each iteration, by default None
        """
        assert self.miller_set is not None

        # Prepare d_params_d_f_calc
        xrs = self.xrs.deep_copy_scatterers()
        fo = self.miller_set.array(
            flex.complex_double(self.miller_set.size(), 0.01)
        ).as_amplitude_array()
        xrs.shake_sites_in_place(rms_difference=0.1)
        model = f_model(f_obs=fo, xray_structure=xrs)
        target = model.target_functor()(compute_gradients=True)
        dpdfc = target.d_target_d_f_calc_work().data()

        # Run
        r = range(n)
        if show_pbar:
            r = tqdm(
                r,
                desc=self.__class__.__name__.ljust(32),
                position=pbar_position,
                leave=pbar_leave,
            )
        for _ in r:
            self.f_calc()
            self.grads(dpdfc)
            if overall_pbar is not None:
                overall_pbar.update()

    def summary(self, header: bool = True, pad: int = None) -> str:
        """Print a summary of the runs

        Parameters
        ----------
        header : bool, optional
            Whether to print a header with the class name, by default True
        pad : int, optional
            padding for each line, by default None

        Returns
        -------
        str
            multi-line string with statistics
        """
        res: dict[str, list] = {}
        for r in self.runtimes:
            l = res.get(r.desc, [])
            l.append(r.time)
            res[r.desc] = l
        out = []
        if header:
            out.append(self.__class__.__name__.ljust(30))
        if pad is None:
            pad = max(len(k) for k in res.keys()) + 2
        assert isinstance(pad, int)
        for key, times in res.items():
            out.append(
                f"{key.rjust(pad)}: {np.mean(times) :.2e} ± {np.std(times) :.2e}"
            )
        return "\n".join(out)

    def save(self, filepath: Union[Path, str] = None):
        """Save results to csv file by appending. Creates the file with a header if it does not exist.

        Parameters
        ----------
        filepath : Path | str, optional
            Where to save the data.
            If str, uses  `Path(__file__).parent / "data" / filepath`.
            If None, uses `Path(__file__).parent / "data" / "runtime.csv"`, by default None
        """
        if filepath is None:
            filepath = "runtime.csv"
        if isinstance(filepath, str):
            filepath = Path(__file__).parent / "data" / filepath
        if not filepath.exists():
            filepath.parent.mkdir(parents=True, exist_ok=True)
            filepath.write_text(
                ",".join(dataclasses.asdict(self.runtimes[0]).keys()) + "\n"
            )
        with filepath.open("a") as f:
            f.write("\n".join(r.csv() for r in self.runtimes) + "\n")

    def __enter__(self) -> "RuntimeBase":
        """Simple context manager to ensure results are saved"""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Saves the result upon exit from context manager"""
        # Save regardless of keyboard interrupt or such
        if len(self.runtimes) > 0:
            self.save()


def cctbx_factory(
    name: str,
    algorithm: str,
    cos_sin_table: bool = False,
) -> type[RuntimeBase]:
    """Generate a RuntimeBase subclass for running cctbx

    Parameters
    ----------
    name : str
        Class name
    algorithm : str
        Sent to cctbx, either of "direct" or "fft"
    cos_sin_table : bool, optional
        Sent to cctbx, by default False

    Returns
    -------
    type[RuntimeBase]
        New RuntimeBase subclass with __name__ attribute set
    """

    class CctbxRunner(RuntimeBase):
        @measure
        def f_calc(self):
            structure_factors.from_scatterers(
                self.miller_set,
                cos_sin_table=cos_sin_table,
            )(
                self.xrs,
                self.miller_set,
                algorithm,
            )

        @measure
        def grads(self, d_target_d_f_calc):
            structure_factors.gradients(
                self.miller_set,
                cos_sin_table=cos_sin_table,
            )(
                self.xrs,
                None,
                self.miller_set,
                d_target_d_f_calc,
                self.xrs.n_parameters(),
                algorithm,
            )

    CctbxRunner.__name__ = name

    return CctbxRunner


def discamb_factory(name: str, method: FCalcMethod, **kwargs) -> type[RuntimeBase]:
    """Generate a RuntimeBase subclass to run DiSCaMB

    Parameters
    ----------
    name : str
        Class name
    method : FCalcMethod
        Passed to pydiscamb.DiscambWrapper
    **kwargs:
        Passed to pydiscamb.DiscambWrapper

    Returns
    -------
    type[RuntimeBase]
        New RuntimeBase subclass with __name__ attribute set

    """

    class DiscambRunner(RuntimeBase):
        @measure
        def f_calc(self):
            w = DiscambWrapper(self.xrs, method, **kwargs)
            w.f_calc(self.miller_set)

        @measure
        def grads_setup(self) -> DiscambWrapper:
            w = DiscambWrapper(self.xrs, method, **kwargs)
            w.set_indices(self.miller_set.indices())
            return w

        @measure
        def grads(self, d_target_d_f_calc):
            w = self.grads_setup()
            w.d_target_d_params(list(d_target_d_f_calc))

    DiscambRunner.__name__ = name
    return DiscambRunner


def get_dispersed_tyrosines(n_tyrosines: int) -> structure:
    """Place `n_tyrosines` tyrosines randomly around in a big P1 box.
    Each tyrosine has a random orientation, and should be far enough away from others to not cause confusion with e.g. bonds or atom typing.

    Parameters
    ----------
    n_tyrosines : int
        Number of tyrosines to disperse

    Returns
    -------
    structure
        Big P1 structure with tyrosines
    """
    trsn = tyrosine().cubic_unit_cell_around_centered_scatterers(3)
    a, _, _, _, _, _ = trsn.unit_cell().parameters()

    # Use default random structure to get sites.
    # Replace each site with tyrosine
    xrs = cctbx_random_structure.xray_structure(
        space_group_symbol="P1",
        n_scatterers=n_tyrosines,
        min_distance=1,
        elements="const",
        general_positions_only=False,
    )
    shifts = xrs.sites_cart() * a
    xrs = reduce(lambda a, b: a.concatenate(b), [trsn for _ in range(n_tyrosines)])

    trsn_sites = trsn.sites_cart()
    sites = flex.vec3_double(0)
    mt = flex.mersenne_twister()  # To rotate the structures
    for shift in shifts:
        # Add random rotation as well
        rot = mt.random_double_r3_rotation_matrix()
        trsn.set_sites_cart(rot * trsn_sites + shift)
        sites = sites.concatenate(trsn.sites_cart())
    xrs.set_sites_cart(sites)
    xrs = xrs.orthorhombic_unit_cell_around_centered_scatterers(1)
    xrs.scattering_type_registry(table="electron")
    xrs.make_scatterer_labels_shelx_compatible_in_place()
    return xrs


def tyrosine() -> structure:
    import iotbx.pdb
    import mmtbx.model
    from libtbx.utils import null_out

    pdb_str = """
CRYST1   16.170   14.591   15.187  90.00  90.00  90.00 P 1
SCALE1      0.061843  0.000000  0.000000        0.00000
SCALE2      0.000000  0.068535  0.000000        0.00000
SCALE3      0.000000  0.000000  0.065846        0.00000
ATOM      1  N   TYR A   4       8.357   9.217   8.801  1.00 10.55           N
ATOM      2  CA  TYR A   4       9.150   8.055   9.050  1.00 10.24           C
ATOM      3  C   TYR A   4      10.419   8.399   9.804  1.00  9.86           C
ATOM      4  O   TYR A   4      10.726   9.591  10.050  1.00 11.39           O
ATOM      5  CB  TYR A   4       9.496   7.352   7.737  1.00 30.00           C
ATOM      6  CG  TYR A   4       8.296   6.791   7.006  1.00 30.00           C
ATOM      7  CD1 TYR A   4       7.820   5.517   7.293  1.00 30.00           C
ATOM      8  CD2 TYR A   4       7.642   7.534   6.034  1.00 30.00           C
ATOM      9  CE1 TYR A   4       6.724   5.000   6.628  1.00 30.00           C
ATOM     10  CE2 TYR A   4       6.545   7.024   5.364  1.00 30.00           C
ATOM     11  CZ  TYR A   4       6.091   5.758   5.665  1.00 30.00           C
ATOM     12  OH  TYR A   4       5.000   5.244   5.000  1.00 30.00           O
ATOM     13  OXT TYR A   4      11.170   7.502  10.187  1.00  9.86           O
ATOM     14  H1  TYR A   4       8.254   9.323   7.923  1.00 10.55           H
ATOM     15  H2  TYR A   4       7.560   9.120   9.184  1.00 10.55           H
ATOM     16  H3  TYR A   4       8.764   9.932   9.140  1.00 10.55           H
ATOM     17  HA  TYR A   4       8.637   7.441   9.599  1.00 10.24           H
ATOM     18  HB2 TYR A   4       9.929   7.989   7.148  1.00 30.00           H
ATOM     19  HB3 TYR A   4      10.097   6.615   7.928  1.00 30.00           H
ATOM     20  HD1 TYR A   4       8.245   5.005   7.942  1.00 30.00           H
ATOM     21  HD2 TYR A   4       7.946   8.389   5.830  1.00 30.00           H
ATOM     22  HE1 TYR A   4       6.415   4.146   6.829  1.00 30.00           H
ATOM     23  HE2 TYR A   4       6.116   7.532   4.714  1.00 30.00           H
ATOM     24  HH  TYR A   4       4.712   5.804   4.444  1.00 30.00           H
"""
    pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
    model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
    xrs = model.get_xray_structure()
    xrs.scattering_type_registry(table="electron")
    return xrs


def sanitize_pdb(file: Path) -> Path:
    import subprocess

    stage_1 = Path(file.name[:4] + "_1_no_disorder").with_suffix(".pdb")
    stage_2 = Path(file.name[:4] + "_2_hydrogens").with_suffix(".pdb")
    stage_3 = Path(file.name[:4] + "_3_no_water").with_suffix(".pdb")
    stage_4 = Path(file.name[:4] + "_4_no_ions").with_suffix(".pdb")
    out = Path(file.name[:4]).with_suffix(".pdb")
    if out.exists() and out.read_text():
        return out

    # Remove disorder
    subprocess.call(
        [
            "phenix.pdbtools",
            str(file),
            "remove_alt_confs=True",
            "output.filename=" + str(stage_1),
        ],
    )

    # Add hydrogens
    with stage_2.open("w") as stdout:
        subprocess.call(
            [
                "phenix.reduce",
                "-DROP_HYDROGENS_ON_ATOM_RECORDS",
                str(stage_1),
            ],
            stdout=stdout,
        )

    # Remove waters
    subprocess.call(
        [
            "phenix.ready_set",
            "pdb_file_name=" + str(stage_2),
            "ligands=False",
            "remove_waters=True",
            "output_file_name=" + str(stage_3.with_suffix("")),
        ],
    )

    # Remove ions
    with stage_4.open("w") as fout, stage_3.open("r") as fin:
        for line in fin:
            if line[-2] == "-":
                line = line[:-3] + "\n"
            fout.write(line)

    out.write_text(stage_4.read_text())

    return out


def download_pdb(code: str) -> Path:
    out = Path(code + "_raw").with_suffix(".pdb")
    if out.exists():
        return out
    import requests

    assert isinstance(code, str)
    assert len(code) == 4

    response = requests.get(f"https://files.rcsb.org/view/{code}.pdb")
    if response.status_code == 404:
        raise ValueError("pdb code not found on rcsb.org")
    elif response.status_code != 200:
        raise RuntimeError("Communication error with rcsb.org")
    pdb_str = response.content.decode("utf-8")
    pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
    model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())

    with out.open("w") as f:
        f.write(model.as_pdb_or_mmcif_string("pdb"))
    assert out.exists()
    return out


def get_structure_from_pdb(code: str) -> structure:

    original = download_pdb(code)
    sanitized = sanitize_pdb(original)

    pdb_str = sanitized.read_text()
    pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
    model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
    xrs = model.get_xray_structure()

    # Produce type file
    w = DiscambWrapper(xrs, FCalcMethod.TAAM, assignment_info=code + ".log")
    type_assignment = w.atom_type_assignment

    model.set_occupancies(
        flex.double(
            [
                0 if type == "" else 0.5 if lcs == "" else 1
                for type, lcs in type_assignment.values()
            ]
        )
    )
    with open(code + "_typed.pdb", "w") as f:
        f.write(model.as_pdb_or_mmcif_string("pdb"))

    return xrs


CctbxFFTWithCosSin = cctbx_factory(
    "CctbxFFTWithCosSin",
    "fft",
    True,
)
CctbxFFTNoCosSin = cctbx_factory(
    "CctbxFFTNoCosSin",
    "fft",
    False,
)
CctbxDirectNoCosSin = cctbx_factory(
    "CctbxDirectNoCosSin",
    "direct",
    False,
)
CctbxDirectWithCosSin = cctbx_factory(
    "CctbxDirectWithCosSin",
    "direct",
    True,
)
# TODO multithreaded cctbx exists right?

DiscambIamStandard = discamb_factory(
    "DiscambIamStandard",
    FCalcMethod.IAM,
)
DiscambIamMacromol = discamb_factory(
    "DiscambIamMacromol",
    FCalcMethod.IAM,
    algorithm="macromol",
)
DiscambIamMultithread = discamb_factory(
    "DiscambIamMultithread",
    FCalcMethod.IAM,
    n_cores=20,
)
DiscambTaamStandard = discamb_factory(
    "DiscambTaamStandard",
    FCalcMethod.TAAM,
    algorithm="standard",
)
DiscambTaamMacromol = discamb_factory(
    "DiscambTaamMacromol",
    FCalcMethod.TAAM,
    algorithm="macromol",
)
DiscambTaamStandardMultithread = discamb_factory(
    "DiscambTaamStandardMultithread",
    FCalcMethod.TAAM,
    algorithm="standard",
    n_cores=20,
)
DiscambTaamMacromolMultithread = discamb_factory(
    "DiscambTaamMacromolMultithread",
    FCalcMethod.TAAM,
    algorithm="macromol",
    n_cores=20,
)
# DiscambTaamGpu = discamb_factory(FCalcMethod.TAAM, # TODO


def full_sweep_main():

    # A little awkward, but this generates around 120k reflections.
    # Taken from a lysozyme crystal (7DER on the PDB)
    miller_set = crystal.symmetry(
        unit_cell=(78.777, 78.777, 37.093, 90, 90, 90), space_group=96
    ).build_miller_set(False, d_min=0.8)

    # Prepare miller sets from around 40k to around 120k reflections
    sets = []
    while miller_set.size() > 40_000:
        sel = flex.bool(miller_set.size(), True)
        for i in range(10_000):
            j = i * miller_set.size() / 10_000
            sel[int(j)] = False
        miller_set = miller_set.select(sel)
        sets.append(miller_set)

    # Number of tyrosines to disperse in a large P1 unit cell
    ns = [1, 5, 10, 20, 50, 100, 150, 200, 250, 300, 400, 500, 750]

    # Number of runs for each parameter set
    n_runs = 5

    pbar = tqdm(
        desc="Overall",
        total=len(ns) * len(sets) * n_runs * len(CALCS),
        position=0,
    )
    for miller_set in sets:
        for n_trsn in ns:
            # Get structure of randomly dispersed tyrosines
            xrs = get_dispersed_tyrosines(n_trsn)
            # Ensure we compute gradients
            xrs.scatterers().flags_set_grads(state=True)
            # Make the indices compatible with the structure
            miller_set = miller.set(
                xrs.crystal_symmetry(),
                miller_set.indices(),
                False,
            )
            # Compute values for progress bar
            n = xrs.scatterers().size()
            hkl = miller_set.size()
            type_assignment = DiscambWrapper(xrs, FCalcMethod.TAAM).atom_type_assignment
            typed = sum(map(lambda el: bool(el[1]), type_assignment.values())) / n
            # Run
            for Calc in tqdm(
                [
                    CctbxFFTWithCosSin,
                    CctbxDirectWithCosSin,
                    # DiscambIamMacromol,
                    DiscambTaamMacromol,
                ],
                position=1,
                desc=f"{n = }, {hkl = }, {typed = :.1%}",
                leave=False,
            ):
                with Calc(xrs, miller_set=miller_set) as calc:
                    try:
                        calc.run(
                            n_runs,
                            show_pbar=True,
                            pbar_position=2,
                            pbar_leave=False,
                            overall_pbar=pbar,
                        )
                    finally:
                        calc.save("full_sweep.csv")


def small_sweep_main():
    # Prepare miller sets from around 40k to around 120k reflections
    def _sets():
        i = 4
        while i**3 < 100_000:
            r1 = list(range((-i) // 2, i // 2))
            i += 1
            r2 = list(range((-i) // 2, i // 2))
            yield flex.miller_index([(h, k, l) for h in r1 for k in r1 for l in r1])
            yield flex.miller_index([(h, k, l) for h in r2 for k in r1 for l in r1])
            yield flex.miller_index([(h, k, l) for h in r2 for k in r2 for l in r1])

    sets = list(_sets())

    # Number of tyrosines to disperse in a large P1 unit cell
    ns = range(1, 400, 5)

    static_n = 1
    static_set = sets[18]

    print("Static number of scatterers: ", static_n * 24)
    print("Static number of reflections:", static_set.size())
    print("Size of scatterer sweep: ", len(ns))
    print("Size of reflection sweep:", len(sets))
    print("Max scatterers:", ns[-1] * 24)
    print("Max reflections:", sets[-1].size())

    # Number of runs for each parameter set
    n_runs = 5
    calcs = [
        CctbxFFTNoCosSin,
        CctbxDirectWithCosSin,
        DiscambIamMacromol,
        DiscambTaamMacromol,
    ]

    pbar = tqdm(
        desc="Overall",
        total=(len(ns) + len(sets)) * n_runs * len(calcs),
        position=0,
    )

    xrs = get_dispersed_tyrosines(static_n)
    # Ensure we compute gradients
    xrs.scatterers().flags_set_grads(state=True)
    # Make the indices compatible with the structure
    for miller_set in sets:
        miller_set = miller.set(
            xrs.crystal_symmetry(),
            miller_set,
            True,
        )
        assert miller_set.is_unique_set_under_symmetry()
        for Calc in calcs:
            calc = Calc(xrs, miller_set=miller_set)
            try:
                calc.run(n_runs, overall_pbar=pbar)
            finally:
                calc.save("refl_sweep.csv")

    for n_trsn in ns:
        xrs = get_dispersed_tyrosines(n_trsn)
        # Ensure we compute gradients
        xrs.scatterers().flags_set_grads(state=True)
        # Make the indices compatible with the structure
        miller_set = miller.set(
            xrs.crystal_symmetry(),
            static_set,
            True,
        )
        # Run
        for Calc in calcs:
            calc = Calc(xrs, miller_set=miller_set)
            try:
                calc.run(n_runs, overall_pbar=pbar)
            finally:
                calc.save("sctr_sweep.csv")


def pdb_main():
    entries = [
        # ("4ZNN", 1.41),
        # ("5S5D", 2.90, "xray"),  # Very big
        ("6IPU", 1.99, "xray"),  # Very big and has DNA
        ("3NIR", 0.48, "xray"),  # Crambin
        ("7ETN", 0.82, "xray"),  # Decent small molecule
        ("7DER", 1.03, "xray"),  # Lysozyme
        ("6GER", 2.00, "xray"),  # 8k atoms after adding H
        ("6G1T", 1.90, "xray"),  # dna and protein
        # ("7THB", 1.64, "xray"),  # pretty triple dna vortex
    ]

    # Number of runs for each parameter set
    n_runs = 5

    for pdb_code, d_min, table in entries:
        xrs = get_structure_from_pdb(pdb_code)
        xrs.scattering_type_registry(table=table)
        xrs.scatterers().flags_set_grads(state=True)
        # Compute values for progress bar
        n = xrs.scatterers().size()
        _w = DiscambWrapper(xrs, FCalcMethod.TAAM)
        type_assignment = _w.atom_type_assignment
        _w.set_d_min(d_min)
        hkl = len(_w.hkl)
        typed = sum(map(lambda el: bool(el[1]), type_assignment.values())) / n

        # Run
        for Calc in tqdm(
            [
                CctbxFFTWithCosSin,
                CctbxDirectWithCosSin,
                DiscambIamMacromol,
                DiscambTaamMacromol,
            ],
            position=0,
            desc=f"{pdb_code} | {n = }, {hkl = }, {typed = :.1%}",
            leave=True,
        ):
            calc = Calc(xrs, d_min=d_min)
            try:
                calc.run(n_runs, show_pbar=pdb_code == "6IPU", pbar_position=1)
            finally:
                calc.save("pdb.csv")


if __name__ == "__main__":
    # small_sweep_main()
    pdb_main()
    # full_sweep_main()
