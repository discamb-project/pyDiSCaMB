import pytest
from cctbx.eltbx.chemical_elements import proper_caps_list as elements

import pydiscamb


def test_init(tyrosine):
    w = pydiscamb.DiscambWrapper(tyrosine, pydiscamb.FCalcMethod.TAAM)


def test_f_calc(tyrosine):
    w = pydiscamb.DiscambWrapper(tyrosine, pydiscamb.FCalcMethod.TAAM)
    fc = w.f_calc(5.0)
    assert isinstance(fc[0], complex)


def _get_single_element_structure(element: str):
    from cctbx.development import random_structure as cctbx_random_structure
    from cctbx.sgtbx import space_group_info

    group = space_group_info(19)
    xrs = cctbx_random_structure.xray_structure(
        space_group_info=group,
        elements=[element],
        general_positions_only=False,
        use_u_iso=True,
        random_u_iso=False,
        random_occupancy=False,
    )
    xrs.scattering_type_registry(table="wk1995")
    return xrs


@pytest.mark.parametrize(
    "el",
    elements()[:36],
)
def test_assignment_light_element(el: str):
    xrs = _get_single_element_structure(el)

    w1 = pydiscamb.DiscambWrapper(xrs, pydiscamb.FCalcMethod.IAM)
    fc1 = w1.f_calc(2)
    w2 = pydiscamb.DiscambWrapper(xrs, pydiscamb.FCalcMethod.TAAM)
    fc2 = w2.f_calc(2)
    assert pytest.approx(fc1, rel=0.01) == fc2


@pytest.mark.parametrize(
    "el",
    elements()[36:86],
)
def test_assignment_heavy_element(el: str):
    xrs = _get_single_element_structure(el)

    w1 = pydiscamb.DiscambWrapper(xrs, pydiscamb.FCalcMethod.IAM)
    fc1 = w1.f_calc(2)
    w2 = pydiscamb.DiscambWrapper(xrs, pydiscamb.FCalcMethod.TAAM)
    fc2 = w2.f_calc(2)
    assert pytest.approx(fc1, rel=0.01) == fc2


@pytest.mark.parametrize(
    ["table", "expected_R"], [("electron", 0.15), ("wk1995", 0.04)]
)
def test_f_calc_approx_IAM(tyrosine, table, expected_R):
    tyrosine.scattering_type_registry(table=table)
    iam = pydiscamb.DiscambWrapper(tyrosine, pydiscamb.FCalcMethod.IAM)
    taam = pydiscamb.DiscambWrapper(tyrosine, pydiscamb.FCalcMethod.TAAM)

    d_min = 2

    fc_iam = iam.f_calc(d_min)
    fc_taam = taam.f_calc(d_min)

    R = sum(abs(abs(a) - abs(b)) for a, b in zip(fc_iam, fc_taam)) / sum(
        abs(a) for a in fc_taam
    )
    assert R < expected_R


def test_from_parameters(tyrosine):
    wrapper = pydiscamb.DiscambWrapper(
        tyrosine,
        model="taam",
        unit_cell_charge=1000,
        scale=True,
    )
    wrapper.set_d_min(2)
    fc = wrapper.f_calc()


def test_logging(tyrosine, tmp_path):
    wrapper = pydiscamb.DiscambWrapper(
        tyrosine,
        model="matts",
        assignment_info=str(tmp_path / "assignment.txt"),
        parameters_info=str(tmp_path / "parameters.log"),
        multipole_cif=str(tmp_path / "structure.cif"),
        unit_cell_charge=0,
        scale=False,
    )
    assert (tmp_path / "assignment.txt").exists()
    assert (tmp_path / "parameters.log").exists()
    assert (tmp_path / "structure.cif").exists()


def test_unit_cell_charge_scaling(tyrosine):
    w1 = pydiscamb.DiscambWrapper(
        tyrosine,
        model="matts",
        unit_cell_charge=-1000,
        scale=True,
    )
    w2 = pydiscamb.DiscambWrapper(
        tyrosine,
        model="matts",
        unit_cell_charge=1000,
        scale=True,
    )
    assert not pytest.approx(w1.f_calc(2)) == w2.f_calc(2)


def test_unit_cell_charge_scaling_off(tyrosine):
    w1 = pydiscamb.DiscambWrapper(
        tyrosine,
        model="matts",
        unit_cell_charge=1000,
        scale=False,
    )
    w2 = pydiscamb.DiscambWrapper(
        tyrosine,
        model="matts",
        unit_cell_charge=0,
        scale=False,
    )
    assert pytest.approx(w1.f_calc(2)) == w2.f_calc(2)


def test_invalid_bank(tyrosine):
    with pytest.raises(
        RuntimeError,
        match="Problem with accessing/openning file 'non-existent bank file' for reading UBDB type bank",
    ):
        w = pydiscamb.DiscambWrapper(
            tyrosine,
            model="matts",
            bank_path="non-existent bank file",
        )


def test_switching_banks(tyrosine):
    for bank in pydiscamb.get_TAAM_databanks():
        w = pydiscamb.DiscambWrapper(
            tyrosine,
            model="matts",
            bank_path=bank,
        )
