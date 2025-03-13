#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/UnitCell.h"
#include "discamb/CrystalStructure/SpaceGroup.h"
#include "discamb/CrystalStructure/SpaceGroupOperation.h"

#include "discamb/BasicUtilities/discamb_version.h"

#include "PythonInterface.hpp"
#include "scattering_table.hpp"
#include "tests.hpp"
#include "assert.hpp"


namespace py = pybind11;

using namespace std;
using namespace discamb;


PYBIND11_MODULE(_cpp_module, m) {
    m.doc() = R"pbdoc(
        DiSCaMB wrapper using pybind11
        -----------------------

        .. currentmodule:: pydiscamb

        .. autosummary::
           :toctree: _generate

    )pbdoc";

    py::register_exception<AssertionError>(m, "DiscambAssertionError", PyExc_AssertionError);

    m.def(
        "get_discamb_version", 
        &discamb_version::version, 
        R"pbdoc(Get the version string for DiSCaMB)pbdoc"
    );

    py::class_<FCalcDerivatives>(m, "FCalcDerivatives")
        .def_readwrite("hkl", &FCalcDerivatives::hkl)
        .def_readwrite("structure_factor", &FCalcDerivatives::structure_factor)
        .def_property_readonly("site_derivatives", &FCalcDerivatives::siteDerivatives)
        .def_readwrite("adp_derivatives", &FCalcDerivatives::adpDerivatives)
        .def_readwrite("occupancy_derivatives", &FCalcDerivatives::occupancyDerivatives)
    ;

    py::class_<TargetFunctionAtomicParamDerivatives>(m, "TargetDerivatives")
        .def_property_readonly("site_derivatives", [](const TargetFunctionAtomicParamDerivatives &t) {
            return vector<double> {
                t.atomic_position_derivatives.x, 
                t.atomic_position_derivatives.y, 
                t.atomic_position_derivatives.z
            };
        })
        .def_readwrite("adp_derivatives", &TargetFunctionAtomicParamDerivatives::adp_derivatives)
        .def_readwrite("occupancy_derivatives", &TargetFunctionAtomicParamDerivatives::occupancy_derivatives)
    ;

    py::class_<PythonInterface>(m, 
            "PythonInterface", 
            R"pbdoc(Calculate structure factors using DiSCaMB)pbdoc"
        )
        .def(py::init<py::object, py::dict>(), py::arg("structure"), py::arg("kwargs"))
        .def(
            "f_calc", 
            &PythonInterface::f_calc, 
            R"pbdoc(Calculate the structure factors for previously set hkl)pbdoc"
        )
        .def(
            "d_f_calc_d_params",
            &PythonInterface::d_f_calc_d_params,
            R"pbdoc(Calculate the structure factors, and derivatives, for previously set hkl)pbdoc"
        )
        .def(
            "d_f_calc_hkl_d_params",
            &PythonInterface::d_f_calc_hkl_d_params,
            R"pbdoc(Calculate the structure factors, and derivatives, for a given hkl)pbdoc",
            py::arg("h"),
            py::arg("k"),
            py::arg("l")
        )
        .def(
            "d_target_d_params",
            &PythonInterface::d_target_d_params,
            py::return_value_policy::take_ownership,
            R"pbdoc(Calculate the derivatives of a target function)pbdoc",
            py::arg("d_target_d_f_calc")
        )
        .def(
            "set_indices",
            &PythonInterface::set_indices,
            R"pbdoc(Set indices for calculating f_calc. Input must be iterable of tuples with three ints)pbdoc",
            py::arg("indices")
        )
        .def(
            "set_d_min",
            &PythonInterface::set_d_min,
            R"pbdoc(Set minimum d-spacing for calculating f_calc)pbdoc",
            py::arg("d_min")
        )
        .def(
            "update_structure",
            &PythonInterface::update_structure,
            R"pbdoc(Update atoms read from the given structure)pbdoc",
            py::arg("structure")
        )
        .def_property_readonly(
            "hkl",
            [](const PythonInterface &i)
            {
                vector<py::tuple> out;
                for (auto el : i.hkl) out.push_back(py::make_tuple(el.x, el.y, el.z));
                return out;
            },
            R"pbdoc(indices)pbdoc"
        )
    ;

    py::class_<GaussianScatteringParameters>(m, "GaussianScatteringParameters")
        .def_readwrite("a", &GaussianScatteringParameters::a)
        .def_readwrite("b", &GaussianScatteringParameters::b)
        .def_readwrite("c", &GaussianScatteringParameters::c)
        .def("__repr__", &GaussianScatteringParameters::repr)
    ;

    m.def(
        "get_table",
        &get_table,
        R"pbdoc(Get a dict of all scatterers in a given table)pbdoc",
        py::arg("table")
    );

    m.def(
        "table_alias",
        &table_alias,
        R"pbdoc(Get DiSCaMB's name for a scattering table. If not found, the input is returned)pbdoc",
        py::arg("table")
    );

    auto m_tests = m.def_submodule("wrapper_tests", "Tests for the wrapper, written in C++");
    m_tests.def(
        "f_calc_custom_gaussian_parameters", 
        &f_calc_custom_gaussian_parameters,
        R"pbdoc(Calculate structure factors for a given structure and reflections, using explicitly provided gaussian scattering parameters)pbdoc",
        py::arg("structure"),
        py::arg("hkl"),
        py::arg("atom_labels"),
        py::arg("a"),
        py::arg("b"),
        py::arg("c")
    );
}
