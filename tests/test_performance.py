import pytest


@pytest.mark.slow
def test_get_crystal_performance(large_random_structure):
    from time import perf_counter
    from taam_sf import DiscambWrapperTests

    n_iter = 1_000

    wrapper = DiscambWrapperTests(large_random_structure)

    start = perf_counter()
    wrapper.test_get_crystal(n_iter)
    end = perf_counter()

    print(
        f"Runtime for {n_iter = :3_}, {len(large_random_structure.scatterers())} scatterers: {end - start :.1f}s"
    )


@pytest.mark.slow
def test_update_atoms_performance(large_random_structure):
    from time import perf_counter
    from taam_sf import DiscambWrapperTests

    n_iter = 10_000

    wrapper = DiscambWrapperTests(large_random_structure)

    start = perf_counter()
    wrapper.test_update_atoms(n_iter)
    end = perf_counter()

    print(
        f"Runtime for {n_iter = :3_}, {len(large_random_structure.scatterers())} scatterers: {end - start :.1f}s"
    )


@pytest.mark.slow
def test_f_calc_IAM_performance(large_random_structure):
    from time import perf_counter
    from taam_sf import DiscambWrapper

    n_iter = 100
    d_min = 4.0

    wrapper = DiscambWrapper(large_random_structure)

    start = perf_counter()
    for _ in range(n_iter):
        sf = wrapper.f_calc_IAM(d_min)
    end = perf_counter()

    print(
        f"Runtime for {n_iter = :3_}, {len(large_random_structure.scatterers())} scatterers, {len(sf)} reflections: {end - start :.1f}s"
    )


@pytest.mark.slow
def test_f_calc_TAAM_performance(lysosyme):
    from time import perf_counter
    from taam_sf import DiscambWrapper

    n_iter = 100
    d_min = 3.0

    wrapper = DiscambWrapper(lysosyme)

    start = perf_counter()
    for _ in range(n_iter):
        sf = wrapper.f_calc_TAAM(d_min)
    end = perf_counter()

    print(
        f"Runtime for {n_iter = :3_}, {len(lysosyme.scatterers())} scatterers, {len(sf)} reflections: {end - start :.1f}s"
    )


@pytest.mark.slow
def test_f_calc_TAAM_interactive_performance(lysosyme):
    from time import perf_counter
    from taam_sf import InteractiveDiscambWrapper, FCalcMethod

    n_iter = 100
    d_min = 3.0

    wrapper = InteractiveDiscambWrapper(lysosyme, d_min, FCalcMethod.TAAM)

    start = perf_counter()
    for _ in range(n_iter):
        sf = wrapper.f_calc()
    end = perf_counter()

    print(
        f"Runtime for {n_iter = :3_}, {len(lysosyme.scatterers())} scatterers, {len(sf)} reflections: {end - start :.1f}s"
    )


@pytest.mark.slow
def test_f_calc_IAM_interactive_performance(large_random_structure):
    from time import perf_counter
    from taam_sf import InteractiveDiscambWrapper, FCalcMethod

    n_iter = 100
    d_min = 4.0

    wrapper = InteractiveDiscambWrapper(large_random_structure, d_min, FCalcMethod.IAM)

    start = perf_counter()
    for _ in range(n_iter):
        sf = wrapper.f_calc()
    end = perf_counter()

    print(
        f"Runtime for {n_iter = :3_}, {len(large_random_structure.scatterers())} scatterers, {len(sf)} reflections: {end - start :.1f}s"
    )
