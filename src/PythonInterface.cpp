#include "PythonInterface.hpp"

#include "discamb/MathUtilities/Vector3.h"
#include "read_structure.hpp"

using namespace std;
using namespace discamb;

namespace py = pybind11;

PythonInterface::PythonInterface(py::object &structure, py::dict kwargs)
    : DiscambStructureFactorCalculator(
          kwargs.cast<nlohmann::json>(), crystal_from_xray_structure(structure),
          anomalous_from_xray_structure(structure)) {};

void PythonInterface::set_indices(py::object &indices) {
    hkl.clear();
    for (auto hkl_py_auto : indices) {
        py::tuple hkl_py = hkl_py_auto.cast<py::tuple>();
        hkl.push_back(Vector3i{hkl_py[0].cast<int>(),
                               hkl_py[1].cast<int>(),
                               hkl_py[2].cast<int>()});
    }
}

void PythonInterface::update_structure(py::object &structure) {
    update_crystal_from_xray_structure(crystal, structure);
    update_anomalous_from_xray_structure(anomalous, structure);
    update_calculator();
}
