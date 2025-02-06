import pytest

import pydiscamb


def test_init(tyrosine):
    w = pydiscamb.DiscambWrapper(tyrosine, pydiscamb.FCalcMethod.TAAM)


def test_f_calc(tyrosine):
    w = pydiscamb.DiscambWrapper(tyrosine, pydiscamb.FCalcMethod.TAAM)
    fc = w.f_calc(5.0)
    assert isinstance(fc, list)


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
    wrapper = pydiscamb.DiscambWrapper.from_TAAM_parameters(
        tyrosine,
        False,
        pydiscamb.taam_parameters.get_default_databank(),
        "",  # No logging
        "",  # No logging
        "",  # No logging
        1000,
        True,
    )
    wrapper.set_d_min(2)
    fc = wrapper.f_calc()


def test_logging(tyrosine, tmp_path):
    wrapper = pydiscamb.DiscambWrapper.from_TAAM_parameters(
        tyrosine,
        False,
        pydiscamb.taam_parameters.get_default_databank(),
        str(tmp_path / "assignment.txt"),
        str(tmp_path / "parameters.log"),
        str(tmp_path / "structure.cif"),
        0,
        False,
    )
    assert (tmp_path / "assignment.txt").exists()
    assert (tmp_path / "parameters.log").exists()
    assert (tmp_path / "structure.cif").exists()


def test_unit_cell_charge_scaling(tyrosine):
    bank = pydiscamb.taam_parameters.get_default_databank()
    w1 = pydiscamb.DiscambWrapper.from_TAAM_parameters(
        tyrosine, False, bank, "", "", "", -1000, True
    )
    w2 = pydiscamb.DiscambWrapper.from_TAAM_parameters(
        tyrosine, False, bank, "", "", "", 1000, True
    )
    assert not pytest.approx(w1.f_calc(2)) == w2.f_calc(2)


def test_unit_cell_charge_scaling_off(tyrosine):
    bank = pydiscamb.taam_parameters.get_default_databank()
    w1 = pydiscamb.DiscambWrapper.from_TAAM_parameters(
        tyrosine, False, bank, "", "", "", 1000, False
    )
    w2 = pydiscamb.DiscambWrapper.from_TAAM_parameters(
        tyrosine, False, bank, "", "", "", 0, False
    )
    assert pytest.approx(w1.f_calc(2)) == w2.f_calc(2)
