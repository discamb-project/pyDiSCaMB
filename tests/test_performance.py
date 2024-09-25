import pytest

@pytest.mark.slow
def test_throughput(large_random_structure):
    from time import perf_counter
    import random
    from taam_sf import DiscambWrapper

    n_iter = 1_000

    wrapper = DiscambWrapper(large_random_structure)
    runtime = 0
    for _ in range(n_iter):
        for scatterer in large_random_structure.scatterers():
            scatterer.site = (random.random(), random.random(), random.random())
        start = perf_counter()
        wrapper._test_get_crystal()
        end = perf_counter()
        runtime += end - start

    print(f"Runtime for {n_iter = :3_}, {len(large_random_structure.scatterers())} scatterers: {runtime :.1f}s")
