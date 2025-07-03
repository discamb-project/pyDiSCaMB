class DiscambAssertionError(AssertionError): ...

class FCalcDerivatives:
    adp_derivatives: list[list[complex]]
    hkl: list[int]
    occupancy_derivatives: list[complex]
    structure_factor: complex
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    @property
    def site_derivatives(self) -> list[list[complex]]: ...

class GaussianScatteringParameters:
    a: list[float]
    b: list[float]
    c: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class PythonInterface:
    def __init__(self, structure: object, kwargs: dict) -> None:
        """__init__(self: pydiscamb._cpp_module.PythonInterface, structure: object, kwargs: dict) -> None"""
    def d_f_calc_d_params(self) -> list[FCalcDerivatives]:
        """d_f_calc_d_params(self: pydiscamb._cpp_module.PythonInterface) -> list[pydiscamb._cpp_module.FCalcDerivatives]

        Calculate the structure factors, and derivatives, for previously set hkl
        """
    def d_f_calc_hkl_d_params(self, h: int, k: int, l: int) -> FCalcDerivatives:
        """d_f_calc_hkl_d_params(self: pydiscamb._cpp_module.PythonInterface, h: int, k: int, l: int) -> pydiscamb._cpp_module.FCalcDerivatives

        Calculate the structure factors, and derivatives, for a given hkl
        """
    def d_target_d_params(self, d_target_d_f_calc: list[complex]) -> list[TargetDerivatives]:
        """d_target_d_params(self: pydiscamb._cpp_module.PythonInterface, d_target_d_f_calc: list[complex]) -> list[pydiscamb._cpp_module.TargetDerivatives]

        Calculate the derivatives of a target function
        """
    def f_calc(self) -> list[complex]:
        """f_calc(self: pydiscamb._cpp_module.PythonInterface) -> list[complex]

        Calculate the structure factors for previously set hkl
        """
    def selected_d_target_d_params(self, d_target_d_f_calc: list[complex], site: bool, adp: bool, occupancy: bool, fp: bool) -> list[TargetDerivatives]:
        """selected_d_target_d_params(self: pydiscamb._cpp_module.PythonInterface, d_target_d_f_calc: list[complex], site: bool, adp: bool, occupancy: bool, fp: bool) -> list[pydiscamb._cpp_module.TargetDerivatives]

        Calculate the derivatives of a target function
        """
    def set_indices(self, indices: object) -> None:
        """set_indices(self: pydiscamb._cpp_module.PythonInterface, indices: object) -> None

        Set indices for calculating f_calc. Input must be iterable of tuples with three ints
        """
    def update_structure(self, structure: object) -> None:
        """update_structure(self: pydiscamb._cpp_module.PythonInterface, structure: object) -> None

        Update atoms read from the given structure
        """
    @property
    def hkl(self) -> list[tuple]: ...

class TargetDerivatives:
    adp_derivatives: list[float]
    occupancy_derivatives: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    @property
    def site_derivatives(self) -> list[float]: ...

def get_discamb_version() -> str:
    """get_discamb_version() -> str

    Get the version string for DiSCaMB
    """
def get_table(table: str) -> dict[str, GaussianScatteringParameters]:
    """get_table(table: str) -> dict[str, pydiscamb._cpp_module.GaussianScatteringParameters]

    Get a dict of all scatterers in a given table
    """
def table_alias(table: str) -> str:
    """table_alias(table: str) -> str

    Get DiSCaMB's name for a scattering table. If not found, the input is returned
    """
