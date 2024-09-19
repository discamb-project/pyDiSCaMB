import taam_sf

from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info

import pytest

def get_random_structure(space_group: int) -> random_structure.xray_structure:
    group = space_group_info(space_group)
    xrs = random_structure.xray_structure(
            space_group_info       = group,
            elements               = ["C","O","Cd","Fe","H", "N"] * 10,
            general_positions_only = False,
            use_u_iso              = True,
            random_u_iso           = True)
    xrs.scattering_type_registry(table="electron")
    return xrs

@pytest.mark.parametrize("space_group", range(1, 231))
def test_IAM_correctness_random_crystal(space_group: int):
    d_min = 4
    xrs = get_random_structure(space_group)
    fcalc_cctbx = xrs.structure_factors(algorithm="direct", d_min=d_min).f_calc().data()
    fcalc_discamb = taam_sf.test_IAM(xrs, d_min)
    assert [sf.real for sf in fcalc_discamb] == pytest.approx(
        [sf.real for sf in fcalc_cctbx],
        rel=0.3,
    )
    assert [sf.imag for sf in fcalc_discamb] == pytest.approx(
        [sf.imag for sf in fcalc_cctbx],
        rel=0.3,
    )
