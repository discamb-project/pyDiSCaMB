import taam_sf

from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info

import pytest


@pytest.mark.parametrize("space_group", range(1, 231))
@pytest.mark.parametrize("with_adps", ["random_adps", None])
@pytest.mark.parametrize("with_occupancy", ["random_occupancy", None])
def test_IAM_correctness_random_crystal(
    space_group: int, with_adps: bool, with_occupancy: bool
):
    # Use bool() on the params instead of true/false to make pytest output more readable
    d_min = 4
    group = space_group_info(space_group)
    xrs = random_structure.xray_structure(
        space_group_info=group,
        elements=["C", "O", "Cd", "Fe", "N", "H"] * 10,
        general_positions_only=False,
        use_u_iso=True,
        random_u_iso=bool(with_adps),
        random_occupancy=bool(with_occupancy),
    )
    xrs.scattering_type_registry(table="electron")
    fcalc_cctbx = xrs.structure_factors(algorithm="direct", d_min=d_min).f_calc().data()
    fcalc_discamb = taam_sf.test_IAM(xrs, d_min)
    assert fcalc_discamb == pytest.approx(fcalc_cctbx, rel=0.3)
