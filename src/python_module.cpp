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
#include "ManagedDiscambWrapper.hpp"


namespace py = pybind11;

using namespace std;
using namespace discamb;



vector<complex<double>> calculate_structure_factors_TAAM(py::object &structure, double d){
    DiscambWrapper w {structure};
    vector< complex<double> > structureFactors = w.f_calc(d, FCalcMethod::TAAM);
    return structureFactors;
}

vector<complex<double>> calculate_structure_factors_IAM(py::object &structure, double d){
    DiscambWrapper w {structure};
    vector< complex<double> > structureFactors = w.f_calc(d, FCalcMethod::IAM);
    return structureFactors;
}


PYBIND11_MODULE(_taam_sf, m) {
    m.doc() = R"pbdoc(
        DiSCaMB wrapper using pybind11
        -----------------------

        .. currentmodule:: taam_sf

        .. autosummary::
           :toctree: _generate

    )pbdoc";

    m.def(
        "get_discamb_version", 
        &discamb_version::version, 
        R"pbdoc(Get the version string for DiSCaMB)pbdoc"
    );

    m.def(
        "calculate_structure_factors_TAAM", 
        &calculate_structure_factors_TAAM, 
        R"pbdoc(Calculate structure factors for a given structure up to a given d-spacing, using the Independent Atom Model)pbdoc",
        py::arg("structure"),
        py::arg("d_min")
    );
    m.def(
        "calculate_structure_factors_IAM", 
        &calculate_structure_factors_IAM, 
        R"pbdoc(Calculate structure factors for a given structure up to a given d-spacing, using the Transferable Aspherical Atom Model)pbdoc",
        py::arg("structure"),
        py::arg("d_min")
    );

    py::enum_<FCalcMethod>(m, 
            "FCalcMethod", 
            R"pbdoc(Enum for specifying the model for atomic form factor calculations)pbdoc"
        )
        .value("IAM", FCalcMethod::IAM, R"pbdoc(Independent Atom Model)pbdoc")
        .value("TAAM", FCalcMethod::TAAM, R"pbdoc(Transferable Aspherical Atom Model)pbdoc")
        .export_values();

    py::class_<DiscambWrapper>(m, 
            "DiscambWrapper", 
            R"pbdoc(Calculate structure factors using DiSCaMB)pbdoc"
        )
        .def(py::init<py::object>(), py::arg("structure"))
        .def(
            "f_calc_IAM", 
            &DiscambWrapper::f_calc_IAM, 
            R"pbdoc(Calculate the structure factors up to a given d-spacing, using the Independent Atom Model)pbdoc",
            py::arg("d_min")
        )
        .def(
            "f_calc_TAAM", 
            &DiscambWrapper::f_calc_TAAM, 
            R"pbdoc(Calculate the structure factors up to a given d-spacing, using the Transferable Aspherical Atom Model)pbdoc",
            py::arg("d_min")
        )
    ;

    py::class_<DiscambWrapperTests, DiscambWrapper>(m, 
        "DiscambWrapperTests", 
        R"pbdoc(Class for testing the wrapper)pbdoc"
        )
        .def(py::init<py::object>())
        .def("test_get_crystal", &DiscambWrapperTests::test_get_crystal)
        .def("test_update_atoms", &DiscambWrapperTests::test_update_atoms)
    ;

    py::class_<ManagedDiscambWrapper, DiscambWrapper>(m, 
        "ManagedDiscambWrapper", 
        R"pbdoc(Calculate structure factors more efficiently by avoiding unnecessary data transfer and re-computation.)pbdoc"
        )
        .def(
            py::init<py::object, double, FCalcMethod>(), 
            py::arg("structure"), 
            py::arg("d_min"), 
            py::arg("method") = FCalcMethod::IAM
        )
        .def(
            "f_calc", 
            &ManagedDiscambWrapper::f_calc, 
            R"pbdoc(Calculate the structure factors with the parameters supplied to the constructor)pbdoc"
        )
    ;
}
