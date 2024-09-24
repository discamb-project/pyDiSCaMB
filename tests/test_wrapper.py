def test_import():
    import taam_sf

    assert isinstance(taam_sf.get_discamb_version(), str)

def test_init(random_structure):
    from taam_sf import DiscambWrapper
    DiscambWrapper(random_structure)

def test_fcalc(random_structure):
    from taam_sf import DiscambWrapper
    DiscambWrapper(random_structure).f_calc_IAM(4.0)

def test_update_structure(random_structure, lysosyme):
    from taam_sf import DiscambWrapper
    wrapper = DiscambWrapper(random_structure)
    wrapper.update_structure(lysosyme)
