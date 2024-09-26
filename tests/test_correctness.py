import taam_sf

from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info

import pytest


@pytest.mark.slow
@pytest.mark.parametrize("space_group", range(1, 231))
@pytest.mark.parametrize("with_adps", ["random adps", None])
@pytest.mark.parametrize("with_occupancy", ["random occupancy", None])
@pytest.mark.parametrize(
    "atoms",
    ["single weak", "single strong", "many weak", "many strong", "mixed strength"],
)
def test_IAM_correctness_random_crystal(
    space_group: int,
    with_adps: bool,
    with_occupancy: bool,
    atoms: str,
):
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
    d_min = 4
    xrs = random_structure.xray_structure(
        space_group_info=group,
        elements=elements,
        general_positions_only=False,
        use_u_iso=True,
        random_u_iso=bool(with_adps),
        random_occupancy=bool(with_occupancy),
    )
    xrs.scattering_type_registry(table="electron")
    fcalc_cctbx = xrs.structure_factors(algorithm="direct", d_min=d_min).f_calc().data()
    fcalc_discamb = taam_sf.calculate_structure_factors_IAM(xrs, d_min)
    # let pytest handle comparing complex numbers
    assert fcalc_discamb == pytest.approx(fcalc_cctbx, abs=0.3, rel=0.3)


if __name__ == "__main__":
    print("Testing a random crystal with a single weak scatterer in space group 19")
    test_IAM_correctness_random_crystal(19, False, False, "single weak")
    print("Test complete, exiting...")
