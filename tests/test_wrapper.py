import pytest


def test_import():
    import taam_sf

    assert isinstance(taam_sf.get_discamb_version(), str)


def test_init(random_structure):
    from taam_sf import DiscambWrapper

    w = DiscambWrapper(random_structure)


def test_fcalc(random_structure):
    from taam_sf import DiscambWrapper

    w = DiscambWrapper(random_structure)
    sf = w.f_calc_IAM(4.0)
    assert isinstance(sf[0], complex)


def test_update_structure(random_structure):
    from taam_sf import DiscambWrapper

    wrapper = DiscambWrapper(random_structure)
    site = random_structure.scatterers()[0].site
    random_structure.scatterers()[0].site = (site[2], site[1], site[0])


def test_update_structure_recalculate_fcalc(random_structure):
    from taam_sf import DiscambWrapper

    wrapper = DiscambWrapper(random_structure)
    sf_before = wrapper.f_calc_IAM(5)
    site = random_structure.scatterers()[0].site
    random_structure.scatterers()[0].site = (site[2], site[1], site[0])
    sf_after = wrapper.f_calc_IAM(5)
    assert not pytest.approx(sf_before) == sf_after


def test_update_structure_interactive(random_structure):
    from taam_sf import InteractiveDiscambWrapper

    wrapper = InteractiveDiscambWrapper(random_structure, 2.0, "electron-IT")
    sf_before = wrapper.f_calc()
    assert pytest.approx(sf_before) == wrapper.f_calc()
    site = random_structure.scatterers()[0].site
    random_structure.scatterers()[0].site = (site[2], site[1], site[0])
    sf_after = wrapper.f_calc()
    assert not pytest.approx(sf_before) == sf_after
