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

    print(f"Runtime for {n_iter = :3_}, {len(large_random_structure.scatterers())} scatterers: {end - start :.1f}s")

@pytest.mark.slow
def test_update_atoms_performance(large_random_structure):
    from time import perf_counter
    from taam_sf import DiscambWrapperTests

    n_iter = 10_000

    wrapper = DiscambWrapperTests(large_random_structure)

    start = perf_counter()
    wrapper.test_update_atoms(n_iter)
    end = perf_counter()

    print(f"Runtime for {n_iter = :3_}, {len(large_random_structure.scatterers())} scatterers: {end - start :.1f}s")

