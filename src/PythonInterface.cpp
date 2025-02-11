#include "PythonInterface.hpp"

#include "discamb/MathUtilities/Vector3.h"
#include "pybind11_json/pybind11_json.hpp"

#include "read_structure.hpp"

using namespace std;
using namespace discamb;

namespace py = pybind11;

nlohmann::json get_calculator_params(py::object structure, FCalcMethod method){
    nlohmann::json calculator_params;
    switch (method)
    {
    case FCalcMethod::IAM: {
        calculator_params = {
            {"model", "iam"},
            {"electron scattering", false},
            {"table", table_from_xray_structure(structure)},
        };
        break;
    }
    case FCalcMethod::TAAM: {
        calculator_params = {
            {"model", "matts"},
            {"electron scattering", table_from_xray_structure(structure).find("electron") != string::npos},
            {"bank path", py::module::import("pydiscamb.taam_parameters").attr("get_default_databank")().cast<string>()},
        };
        break;
    }
    default:
        break;
    }
    return calculator_params;
}

PythonInterface::PythonInterface(py::object structure) :
    PythonInterface(structure, get_calculator_params(structure, FCalcMethod::IAM))
    {};

PythonInterface::PythonInterface(py::object structure, FCalcMethod method) :
    PythonInterface(structure, get_calculator_params(structure, method))
    {};

PythonInterface::PythonInterface(py::object structure, py::dict kwargs) :
    PythonInterface(structure, kwargs.cast<nlohmann::json>())
    {};
    
PythonInterface::PythonInterface(py::object structure, nlohmann::json calculator_params) :
    mStructure(structure),
    DiscambStructureFactorCalculator(
        calculator_params,
        crystal_from_xray_structure(structure),
        anomalous_from_xray_structure(structure)
    ) 
    {};

void PythonInterface::set_indices(py::object indices){
    hkl.clear();
    for (auto hkl_py_auto : indices){
        py::tuple hkl_py = hkl_py_auto.cast<py::tuple>();
        hkl.push_back(Vector3i {
            hkl_py[0].cast<int>(),
            hkl_py[1].cast<int>(),
            hkl_py[2].cast<int>()
        });
    }
}

void PythonInterface::set_d_min(const double d_min){
    bool anomalous_flag = mStructure.attr("scatterers")().attr("count_anomalous")().cast<int>() != 0;
    py::object miller_py = mStructure.attr("build_miller_set")(anomalous_flag, d_min);

    set_indices(miller_py.attr("indices")());
}

vector<complex<double>> calculate_structure_factors(py::object structure, double d, FCalcMethod method){
    PythonInterface w {structure, method};
    w.set_d_min(d);
    vector< complex<double> > structureFactors = w.f_calc();
    return structureFactors;
}

vector<complex<double>> calculate_structure_factors_TAAM(py::object structure, double d){
    return calculate_structure_factors(structure, d, FCalcMethod::TAAM);
}

vector<complex<double>> calculate_structure_factors_IAM(py::object structure, double d){
    return calculate_structure_factors(structure, d, FCalcMethod::IAM);
}
