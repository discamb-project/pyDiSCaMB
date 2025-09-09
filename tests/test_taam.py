import sys

import pytest

from pydiscamb import DiscambWrapper, FCalcMethod
from pydiscamb.taam_parameters import get_TAAM_databanks, is_MATTS_installed


def test_init(tyrosine):
    w = DiscambWrapper(tyrosine, FCalcMethod.TAAM)


def test_f_calc(tyrosine):
    w = DiscambWrapper(tyrosine, FCalcMethod.TAAM)
    fc = w.f_calc(5.0)
    assert isinstance(fc[0], complex)


@pytest.mark.parametrize(
    ["table", "expected_R"], [("electron", 0.15), ("wk1995", 0.04)]
)
def test_f_calc_approx_IAM(tyrosine, table, expected_R):
    tyrosine.scattering_type_registry(table=table)
    iam = DiscambWrapper(tyrosine, FCalcMethod.IAM)
    taam = DiscambWrapper(tyrosine, FCalcMethod.TAAM)

    d_min = 2

    fc_iam = iam.f_calc(d_min)
    fc_taam = taam.f_calc(d_min)

    R = sum(abs(abs(a) - abs(b)) for a, b in zip(fc_iam, fc_taam)) / sum(
        abs(a) for a in fc_taam
    )
    assert R < expected_R


def test_from_parameters(tyrosine):
    wrapper = DiscambWrapper(
        tyrosine,
        FCalcMethod.TAAM,
        unit_cell_charge=1000,
        scale=True,
    )
    wrapper.set_d_min(2)
    fc = wrapper.f_calc()


def test_logging(tyrosine, tmp_path):
    wrapper = DiscambWrapper(
        tyrosine,
        FCalcMethod.TAAM,
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
    w1 = DiscambWrapper(
        tyrosine,
        FCalcMethod.TAAM,
        unit_cell_charge=-1000,
        scale=True,
    )
    w2 = DiscambWrapper(
        tyrosine,
        FCalcMethod.TAAM,
        unit_cell_charge=1000,
        scale=True,
    )
    assert not pytest.approx(w1.f_calc(2)) == w2.f_calc(2)


def test_unit_cell_charge_scaling_off(tyrosine):
    w1 = DiscambWrapper(
        tyrosine,
        FCalcMethod.TAAM,
        unit_cell_charge=1000,
        scale=False,
    )
    w2 = DiscambWrapper(
        tyrosine,
        FCalcMethod.TAAM,
        unit_cell_charge=0,
        scale=False,
    )
    assert pytest.approx(w1.f_calc(2)) == w2.f_calc(2)


@pytest.mark.xfail(
    condition=sys.platform.startswith("win"),
    reason="Bug on windows where exception object is unpopulated from discamb",
    raises=AssertionError,  # Thrown by pytest when raises fails
)
def test_invalid_bank(tyrosine):
    with pytest.raises(
        RuntimeError,
        match="Problem with accessing/openning file 'non-existent bank file' for reading UBDB type bank",
    ):
        w = DiscambWrapper(
            tyrosine,
            FCalcMethod.TAAM,
            bank_path="non-existent bank file",
        )


@pytest.mark.skipif(not is_MATTS_installed(), reason="Must have MATTS installed")
def test_switching_banks(tyrosine):
    banks = get_TAAM_databanks()
    assert len(banks) > 1
    for bank in banks:
        w = DiscambWrapper(
            tyrosine,
            FCalcMethod.TAAM,
            bank_path=bank,
        )


@pytest.mark.skipif(not is_MATTS_installed(), reason="Must have MATTS installed")
@pytest.mark.parametrize(
    "algorithm",
    [
        "standard",
        "macromol",
    ],
)
def test_frozen_lcs(tyrosine, algorithm):
    from cctbx.array_family import flex

    xrs = tyrosine.deep_copy_scatterers()

    indices = flex.miller_index(10, (1, 2, 3))
    d_target_d_fcalc = [1 + 1j] * 10

    # Get gradients from default parameters
    default = DiscambWrapper(
        xrs, FCalcMethod.TAAM, frozen_lcs=False, algorithm=algorithm
    )
    # Check assignment
    assert all(lcs != "" for atomtype, lcs in default.atom_type_assignment.values())

    default.set_indices(indices)
    grads = default.d_target_d_params(d_target_d_fcalc)
    default_1st_site = grads[0].site_derivatives

    # Get gradients with frozen_lcs
    frozen = DiscambWrapper(xrs, FCalcMethod.TAAM, frozen_lcs=True, algorithm=algorithm)
    assert all(lcs != "" for atomtype, lcs in frozen.atom_type_assignment.values())

    frozen.set_indices(indices)
    grads = frozen.d_target_d_params(d_target_d_fcalc)
    frozen_1st_site = grads[0].site_derivatives

    # Shake the scatterers and update the wrappers
    xrs.shake_sites_in_place(0.2)
    default.update_structure(xrs)
    frozen.update_structure(xrs)

    # Compute gradients again
    grads = default.d_target_d_params(d_target_d_fcalc)
    default_shaken_1st_site = grads[0].site_derivatives

    grads = frozen.d_target_d_params(d_target_d_fcalc)
    frozen_shaken_1st_site = grads[0].site_derivatives

    assert all(a == b for a, b in zip(default_1st_site, frozen_1st_site))

    # This xfails on macromol, raise different error to seperate from other assertions failing
    if not all(a != b for a, b in zip(default_shaken_1st_site, frozen_shaken_1st_site)):
        raise RuntimeError
