#include "PythonInterface.hpp"

#include "discamb/MathUtilities/Vector3.h"

#include "read_structure.hpp"

using namespace std;
using namespace discamb;

namespace py = pybind11;


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

void PythonInterface::update_structure(py::object structure){
    update_crystal_from_xray_structure(crystal, structure);
    update_anomalous_from_xray_structure(anomalous, structure);
    // DiscambStructureFactorCalculator::update_calculator is called before all relevant operations
    // ensuring the updates here are handled
}
