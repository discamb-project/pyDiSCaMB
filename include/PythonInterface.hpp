#pragma once

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/Scattering/SfCalculator.h"

#include "json.hpp"
#include "pybind11_json/pybind11_json.hpp"

#include <vector>
#include <complex>

#include "DiscambStructureFactorCalculator.hpp"

namespace py = pybind11;


class PythonInterface : public DiscambStructureFactorCalculator{
    public:
        PythonInterface(py::object structure, py::dict kwargs) : PythonInterface(structure, kwargs.cast<nlohmann::json>()) {};
        PythonInterface(py::object structure, nlohmann::json calculator_params);

        void set_indices(py::object indices);
        void set_d_min(const double d_min);
        void update_structure(py::object structure);

    private:
        py::object mStructure;
};
