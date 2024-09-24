#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/UnitCell.h"
#include "discamb/CrystalStructure/SpaceGroup.h"
#include "discamb/CrystalStructure/SpaceGroupOperation.h"

#include "discamb/BasicUtilities/discamb_version.h"

#include "DiscambWrapper.hpp"


namespace py = pybind11;

using namespace std;
using namespace discamb;



vector<complex<double>> test_TAAM(const py::object structure, double d){
    // Crystal crystal;
    // DiscambWrapper::structure_to_crystal(structure, crystal);
    // vector<Vector3i> hkl;
    // DiscambWrapper::structure_to_hkl(structure, d, hkl);

    // structureFactors.reserve(hkl.size());
    // calculateSfTaamMinimal(crystal, hkl, structureFactors);

    vector< complex<double> > structureFactors;
    return structureFactors;
}

vector<complex<double>> test_IAM(const py::object structure, double d){
    // Crystal crystal;
    // DiscambWrapper::structure_to_crystal(structure, crystal);
    // vector<Vector3i> hkl;
    // DiscambWrapper::structure_to_hkl(structure, d, hkl);

    // structureFactors.reserve(hkl.size());
    // calculateSfIamMinimal(crystal, hkl, structureFactors);

    vector< complex<double> > structureFactors;
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
        .def("update_structure", &DiscambWrapper::update_structure)
        .def("f_calc_IAM", [](DiscambWrapper &self, double d_min){ return self.f_calc(d_min, FCalcMethod::IAM); })
        .def("f_calc_TAAM", [](DiscambWrapper &self, double d_min){ return self.f_calc(d_min, FCalcMethod::TAAM); })
    ;
}
