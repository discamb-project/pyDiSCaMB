#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/UnitCell.h"
#include "discamb/CrystalStructure/SpaceGroup.h"
#include "discamb/CrystalStructure/SpaceGroupOperation.h"

#include "discamb/BasicUtilities/discamb_version.h"

#include "DiscambWrapper.hpp"
#include "DiscambWrapperTests.hpp"


namespace py = pybind11;

using namespace std;
using namespace discamb;



vector<complex<double>> test_TAAM(py::object &structure, double d){
    DiscambWrapper w {structure};
    vector< complex<double> > structureFactors = w.f_calc(d, FCalcMethod::TAAM);
    return structureFactors;
}

vector<complex<double>> test_IAM(py::object &structure, double d){
    DiscambWrapper w {structure};
    vector< complex<double> > structureFactors = w.f_calc(d, FCalcMethod::IAM);
    return structureFactors;
}


PYBIND11_MODULE(_taam_sf, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: taam_sf

        .. autosummary::
           :toctree: _generate

           get_discamb_version
    )pbdoc";

    m.def("get_discamb_version", &discamb_version::version, R"pbdoc(
        Get the version string for DiSCaMB
    )pbdoc");

    m.def("test_TAAM", &test_TAAM, "placeholder docstring");
    m.def("test_IAM", &test_IAM, "placeholder docstring");

    py::class_<DiscambWrapper>(m, "DiscambWrapper")
        .def(py::init<py::object>())
        .def("f_calc_IAM", &DiscambWrapper::f_calc_IAM)
        .def("f_calc_TAAM", &DiscambWrapper::f_calc_TAAM)
    ;

    py::class_<DiscambWrapperTests, DiscambWrapper>(m, "DiscambWrapperTests")
        .def(py::init<py::object>())
        .def("test_get_crystal", &DiscambWrapperTests::test_get_crystal)
        .def("test_update_atoms", &DiscambWrapperTests::test_update_atoms)
    ;
}
