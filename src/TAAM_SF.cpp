#include <pybind11/pybind11.h>

#include "discamb_wrapper.hpp"


namespace py = pybind11;

PYBIND11_MODULE(_taam_sf, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: taam_sf

        .. autosummary::
           :toctree: _generate

           get_discamb_version
    )pbdoc";

    m.def("get_discamb_version", &get_discamb_version, R"pbdoc(
        Get the version string for DiSCaMB
    )pbdoc");

}
