from cctbx.array_family import flex
from cctbx import miller, xray
from mmtbx import f_model
import pytest


def test_import():
    import pydiscamb

    assert isinstance(pydiscamb.get_discamb_version(), str)


def test_init(random_structure):
    from pydiscamb import DiscambWrapper

    w = DiscambWrapper(random_structure)


def test_fcalc(random_structure):
    from pydiscamb import DiscambWrapper

    w = DiscambWrapper(random_structure)
    sf = w.f_calc(4.0)
    assert isinstance(sf[0], complex)


def test_fcalc_with_d_min(random_structure):
    from pydiscamb import DiscambWrapper

    w = DiscambWrapper(random_structure)
    w.set_d_min(4.0)
    sf = w.f_calc()
    assert isinstance(sf[0], complex)


def test_fcalc_with_indices(random_structure):
    from pydiscamb import DiscambWrapper

    inds = [(0, 1, 2), (0, 1, 2), (10, 20, 30)]

    w = DiscambWrapper(random_structure)
    w.set_indices(inds)
    sf = w.f_calc()
    assert len(sf) == len(inds)
    assert isinstance(sf[0], complex)
    assert pytest.approx(sf[0]) == sf[1]
    assert pytest.approx(sf[0]) != sf[2]


def test_fcalc_with_no_indices(random_structure):
    from pydiscamb import DiscambWrapper

    w = DiscambWrapper(random_structure)
    sf = w.f_calc()
    assert len(sf) == 0


def test_update_structure(random_structure):
    from pydiscamb import DiscambWrapper

    wrapper = DiscambWrapper(random_structure)
    site = random_structure.scatterers()[0].site
    random_structure.scatterers()[0].site = (site[2], site[1], site[0])


@pytest.mark.xfail(reason="The wrapper is no longer updated along with the structure")
def test_update_structure_recalculate_fcalc(random_structure):
    from pydiscamb import DiscambWrapper

    wrapper = DiscambWrapper(random_structure)
    sf_before = wrapper.f_calc(5)
    site = random_structure.scatterers()[0].site
    random_structure.scatterers()[0].site = (site[2], site[1], site[0])
    sf_after = wrapper.f_calc(5)
    assert not pytest.approx(sf_before) == sf_after


@pytest.mark.parametrize(
    ["p", "dp"],
    [
        (False, True),
        (True, False),
        (True, True),
    ],
)
def test_anomalous_scattering(p: bool, dp: bool, random_structure):
    from pydiscamb import DiscambWrapper

    wrapper = DiscambWrapper(random_structure)
    sf_before = wrapper.f_calc(5)

    sf_before = random_structure.structure_factors(d_min=5).f_calc().data()
    rb = [sf.real for sf in sf_before]
    ib = [sf.imag for sf in sf_before]
    b = [abs(sf) for sf in sf_before]

    if p:
        random_structure.shake_fps()
    if dp:
        random_structure.shake_fdps()

    wrapper = DiscambWrapper(random_structure)
    sf_after = wrapper.f_calc(5)
    sf_after = random_structure.structure_factors(d_min=5).f_calc().data()
    ra = [sf.real for sf in sf_after]
    ia = [sf.imag for sf in sf_after]
    a = [abs(sf) for sf in sf_after]

    assert not pytest.approx(rb) == ra
    assert not pytest.approx(ib) == ia
    assert not pytest.approx(b) == a


def test_f_calc_type(tyrosine):
    from pydiscamb import DiscambWrapper

    fc_c = tyrosine.structure_factors(d_min=2).f_calc()

    w = DiscambWrapper(tyrosine)
    w.set_indices(fc_c.indices())
    assert isinstance(w.f_calc(), flex.complex_double)

    w = DiscambWrapper(tyrosine)
    assert isinstance(w.f_calc(fc_c), miller.array)


@pytest.mark.parametrize("taam", [True, False])
def test_default_init_with_kwargs(tyrosine, taam):
    from pydiscamb import DiscambWrapper, FCalcMethod

    m = FCalcMethod.TAAM if taam else FCalcMethod.IAM
    w1 = DiscambWrapper(tyrosine, m)
    w2 = DiscambWrapper(tyrosine, m, random_unused_kwarg=42)

    fc1 = w1.f_calc(2.0)
    fc2 = w2.f_calc(2.0)
    assert pytest.approx(list(fc1)) == list(fc2)
