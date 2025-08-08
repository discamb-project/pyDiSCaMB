# Script to measure the runtime of the wrapper, Python and C++ side

from measure_fcalc_gradients_runtime import get_structure_from_pdb
from time import perf_counter

from cctbx.xray.structure_factors import gradients


def main():
    xrs = get_structure_from_pdb("7DER")
    fc = xrs.structure_factors(d_min=1.03).f_calc()
    ms = fc.set()

    start = perf_counter()
    gradients(
        ms,
    )(
        xrs,
        None,
        ms,
        fc.data(),
        xrs.n_parameters(),
        "taam",
    )
    end = perf_counter()

    print("total", end - start)


if __name__ == "__main__":
    main()
