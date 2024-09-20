import taam_sf

from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info

import pytest


@pytest.fixture
def structure() -> random_structure.xray_structure:
    group = space_group_info(19)
    xrs = random_structure.xray_structure(
        space_group_info=group,
        elements=["Au", "C"] * 3,
        general_positions_only=False,
        use_u_iso=True,
        random_u_iso=False,
        random_occupancy=False,
    )
    return xrs


@pytest.mark.parametrize("model", [taam_sf.test_IAM, taam_sf.test_TAAM])
class TestSF:

    def test_is_callable(self, model):
        assert callable(model)

    def test_with_random_structure(self, model, structure):
        model(structure, 3)
