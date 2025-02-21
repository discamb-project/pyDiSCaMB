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

def test_assignment_log_urecognized_structure(tmp_path):
    logfile = tmp_path / "assignment.log"

    from cctbx.development import random_structure as cctbx_random_structure
    from cctbx.sgtbx import space_group_info

    group = space_group_info(19)
    xrs = cctbx_random_structure.xray_structure(
        space_group_info=group,
        elements=elements()[:86],
        general_positions_only=False,
        use_u_iso=True,
        random_u_iso=False,
        random_occupancy=False,
    )
    xrs.scattering_type_registry(table="wk1995")

    w = pydiscamb.DiscambWrapper(xrs, pydiscamb.FCalcMethod.TAAM, assignment_info=str(logfile))

    with logfile.open("r") as f:
        assert f.readline() == "Atom type assigned to 0 of 86.\n"
        assert f.readline() == "Atoms with unassigned atom types :\n"
        for i in range(36, 86):
            assert f.readline() == f"{elements()[i]}{i + 1}\n".rjust(8)
        assert f.readline() == "Atoms with spherical atom type :\n"
        for i in range(36):
            assert f.readline() == f"{elements()[i]}{i + 1}\n".rjust(8)

def test_assignment_log_tyrosine(tyrosine, tmp_path):
    logfile = tmp_path / "assignment.log"

    w = pydiscamb.DiscambWrapper(tyrosine, pydiscamb.FCalcMethod.TAAM, assignment_info=str(logfile))
    
    with logfile.open("r") as f:
        assert f.readline() == "Atom type assigned to 24 of 24.\n"
        assert f.readline() == "Atoms with assigned atom types and local coordinate systems:\n"
        assert f.readline() == '   pdb=" N   TYR A   4 "   N401a    Z pdb=" CA  TYR A   4 " X pdb=" H2  TYR A   4 "\n'
        assert f.readline() == '   pdb=" CA  TYR A   4 "    C414    X pdb=" N   TYR A   4 " Y pdb=" C   TYR A   4 "\n'
        assert f.readline() == '   pdb=" C   TYR A   4 "    C301    Z pdb=" CA  TYR A   4 " X pdb=" OXT TYR A   4 "\n'
        assert f.readline() == '   pdb=" O   TYR A   4 "    O101    X pdb=" C   TYR A   4 " Y pdb=" OXT TYR A   4 "\n'
        assert f.readline() == '   pdb=" CB  TYR A   4 "    C405    X pdb=" CG  TYR A   4 " Y pdb=" CA  TYR A   4 "\n'
        assert f.readline() == '   pdb=" CG  TYR A   4 "    C332    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CD2 TYR A   4 "\n'
        assert f.readline() == '   pdb=" CD1 TYR A   4 "    C330    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CE1 TYR A   4 "\n'
        assert f.readline() == '   pdb=" CD2 TYR A   4 "    C330    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CE2 TYR A   4 "\n'
        assert f.readline() == '   pdb=" CE1 TYR A   4 "    C330    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CZ  TYR A   4 "\n'
        assert f.readline() == '   pdb=" CE2 TYR A   4 "    C330    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CZ  TYR A   4 "\n'
        assert f.readline() == '   pdb=" CZ  TYR A   4 "   C341a    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CE2 TYR A   4 "\n'
        assert f.readline() == '   pdb=" OH  TYR A   4 "    O204    X pdb=" CZ  TYR A   4 " Y pdb=" HH  TYR A   4 "\n'
        assert f.readline() == '   pdb=" OXT TYR A   4 "    O101    X pdb=" C   TYR A   4 " Y pdb=" O   TYR A   4 "\n'
        assert f.readline() == '   pdb=" H1  TYR A   4 "    H105    Z pdb=" N   TYR A   4 " X pdb=" CA  TYR A   4 "\n'
        assert f.readline() == '   pdb=" H2  TYR A   4 "    H105    Z pdb=" N   TYR A   4 " X pdb=" CA  TYR A   4 "\n'
        assert f.readline() == '   pdb=" H3  TYR A   4 "    H105    Z pdb=" N   TYR A   4 " X pdb=" CA  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HA  TYR A   4 "    H103    Z pdb=" CA  TYR A   4 " X pdb=" N   TYR A   4 "\n'
        assert f.readline() == '   pdb=" HB2 TYR A   4 "    H102    Z pdb=" CB  TYR A   4 " X pdb=" CG  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HB3 TYR A   4 "    H102    Z pdb=" CB  TYR A   4 " X pdb=" CG  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HD1 TYR A   4 "    H104    Z pdb=" CD1 TYR A   4 " X pdb=" CE1 TYR A   4 "\n'
        assert f.readline() == '   pdb=" HD2 TYR A   4 "    H104    Z pdb=" CD2 TYR A   4 " X pdb=" CE2 TYR A   4 "\n'
        assert f.readline() == '   pdb=" HE1 TYR A   4 "    H104    Z pdb=" CE1 TYR A   4 " X pdb=" CZ  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HE2 TYR A   4 "    H104    Z pdb=" CE2 TYR A   4 " X pdb=" CZ  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HH  TYR A   4 "    H114    Z pdb=" OH  TYR A   4 " X pdb=" CZ  TYR A   4 "\n'


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
