import pytest


@pytest.mark.slow
def test_get_crystal_performance(large_random_structure):
    from time import perf_counter
    from pydiscamb import DiscambWrapperTests

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
    from pydiscamb import DiscambWrapperTests

    n_iter = 10_000

    wrapper = DiscambWrapperTests(large_random_structure)

    start = perf_counter()
    wrapper.test_update_atoms(n_iter)
    end = perf_counter()

    print(
        f"Runtime for {n_iter = :3_}, {len(large_random_structure.scatterers())} scatterers: {end - start :.1f}s"
    )

def _test_f_calc_performance(structure, method, n_iter, d_min):
    from time import perf_counter
    from pydiscamb import DiscambWrapperTests, FCalcMethod
    method = {"IAM": FCalcMethod.IAM, "TAAM": FCalcMethod.TAAM}[method]

    wrapper = DiscambWrapperTests(structure, method)

    start = perf_counter()
    runtime = wrapper.get_f_calc_runtime(n_iter, d_min)
    mid = perf_counter()
    long_runtime = wrapper.get_f_calc_runtime_with_atom_updates(n_iter, d_min)
    end = perf_counter()

    for _ in range(n_iter):
        structure.structure_factors(d_min=d_min, algorithm="direct").f_calc().data()
    cctbx_end = perf_counter()

    print(
        f"""{method.name} runtime for {n_iter = :3_}, {len(structure.scatterers())} scatterers, {d_min = :.1f} ({len(wrapper.f_calc(d_min))} hkls): 
        No updates:   {runtime / 1000 :.1f}s ({(mid - start) * 1000 - runtime :.1f}ms overhead)
        With updates: {long_runtime / 1000 :.1f}s ({(end - mid) * 1000 - runtime :.1f}ms overhead)
        cctbx (IAM):  {cctbx_end - end :.1f}s
        """
    )

@pytest.mark.slow
@pytest.mark.parametrize("method", ["IAM", "TAAM"])
def test_f_calc_performance_large_structure(lysozyme, method):
    _test_f_calc_performance(lysozyme, method, n_iter=10, d_min=1)

@pytest.mark.slow
@pytest.mark.parametrize("method", ["IAM", "TAAM"])
def test_f_calc_performance_small_structure(tyrosine, method):
    _test_f_calc_performance(tyrosine, method, n_iter=10_000, d_min=3)
