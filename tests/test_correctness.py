import pydiscamb

from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
from cctbx.array_family import flex

import pytest


def compare_structure_factors(f_calc_a, f_calc_b):
    f_calc_a = flex.abs(f_calc_a)
    f_calc_b = flex.abs(f_calc_b)
    scale = flex.sum(f_calc_a * f_calc_b) / flex.sum(f_calc_b * f_calc_b)
    num = flex.sum(flex.abs(f_calc_a - scale * f_calc_b))
    den = flex.sum(flex.abs(f_calc_a + scale * f_calc_b))
    diff = flex.abs(f_calc_a - f_calc_b)
    return num / den * 2 * 100.0  # , flex.mean(diff), flex.max(diff)


def get_random_crystal(
    space_group: int,
    with_adps: bool,
    with_occupancy: bool,
    with_anomalous: str,
    atoms: str,
    scattering_table: str,
) -> random_structure.xray_structure:
    # We use bool() on some of the params instead of true/false to make pytest output more readable
    if atoms == "single weak":
        elements = ["C"]
    if atoms == "single strong":
        elements = ["Au"]
    elif atoms == "many weak":
        elements = ["C", "O", "N", "H"] * 10
    elif atoms == "many strong":
        elements = ["Au", "Ag", "Cd", "Fe"] * 10
    elif atoms == "mixed strength":
        elements = ["Au", "Cd", "C", "O", "H"] * 10
    group = space_group_info(space_group)
    xrs = random_structure.xray_structure(
        space_group_info=group,
        elements=elements,
        general_positions_only=False,
        use_u_iso=with_adps == "random u_iso",
        use_u_aniso=with_adps == "random u_aniso",
        random_u_iso=bool(with_adps),
        random_occupancy=bool(with_occupancy),
    )
    if "fprime" in with_anomalous:
        xrs.shake_fps()
    if "fdoubleprime" in with_anomalous:
        xrs.shake_fdps()
    xrs.scattering_type_registry(table=scattering_table)
    return xrs


def get_IAM_correctness_score(xrs: random_structure.xray_structure) -> float:
    d_min = 2
    fcalc_cctbx = xrs.structure_factors(algorithm="direct", d_min=d_min).f_calc().data()
    fcalc_discamb = pydiscamb.calculate_structure_factors_IAM(xrs, d_min)

    fcalc_discamb = flex.complex_double(fcalc_discamb)
    score = compare_structure_factors(fcalc_cctbx, fcalc_discamb)
    return score


@pytest.mark.slow
@pytest.mark.parametrize(
    "scattering_table",
    ["it1992", "wk1995", "electron"],
)
@pytest.mark.parametrize("space_group", range(1, 231))
@pytest.mark.parametrize("with_adps", ["random u_iso", "random u_aniso", None])
@pytest.mark.parametrize("with_occupancy", ["random occupancy", None])
@pytest.mark.parametrize(
    "with_anomalous",
    ["no anomalous", "fprime", "fdoubleprime", "fprime + fdoubleprime"],
)
@pytest.mark.parametrize(
    "atoms",
    ["single weak", "single strong", "many weak", "many strong", "mixed strength"],
)
def test_IAM_correctness_random_crystal(
    space_group: int,
    with_adps: bool,
    with_occupancy: bool,
    with_anomalous: str,
    atoms: str,
    scattering_table: str,
):
    xrs = get_random_crystal(
        space_group, with_adps, with_occupancy, with_anomalous, atoms, scattering_table
    )
    score = get_IAM_correctness_score(xrs)
    # Use 0.05% as threshold
    assert score < 0.0005

if __name__ == "__main__":
    space_group: int = 19
    with_adps: bool = False
    with_occupancy: bool = False
    with_anomalous: str = ""
    atoms: str = "single weak"
    scattering_table: str = "electron"
    xrs = get_random_crystal(
        space_group, with_adps, with_occupancy, with_anomalous, atoms, scattering_table
    )
    score = get_IAM_correctness_score(xrs)
