def test_import():
    import taam_sf

    assert isinstance(taam_sf.get_discamb_version(), str)
