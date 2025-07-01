# Script to measure the runtime of the wrapper specifically, e.g. initialization and cctbx-discamb translation

from pydiscamb._cpp_module._wrapper_tests import TimedInterface
from pydiscamb._cpp_module import PythonInterface
from pydiscamb.taam_parameters import get_default_databank
from measure_fcalc_gradients_runtime import get_dispersed_tyrosines


def main():
    xrs = get_dispersed_tyrosines(100)
    print(xrs)
    kwargs = {
        "model": "taam",
        "bank path": get_default_databank(),
        "algorithm": "macromol",
    }
    hkl = [(h, k, l) for h in range(30) for k in range(30) for l in range(30)]
    print(len(hkl))

    w = TimedInterface(xrs, kwargs)
    w.set_indices(hkl)
    w.f_calc()
    print(w.get_runtimes())


if __name__ == "__main__":
    main()
