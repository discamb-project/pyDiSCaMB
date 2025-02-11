#pragma once

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/Scattering/SfCalculator.h"

#include "json.hpp"

#include <vector>
#include <complex>

#include "DiscambStructureFactorCalculator.hpp"

namespace py = pybind11;

enum FCalcMethod {
    IAM,
    TAAM
};


class PythonInterface : public DiscambStructureFactorCalculator{
    public:
        PythonInterface(py::object structure);
        PythonInterface(py::object structure, FCalcMethod method);
        PythonInterface(py::object structure, py::dict kwargs);
        PythonInterface(py::object structure, nlohmann::json calculator_params);

        void set_indices(py::object indices);
        void set_d_min(const double d_min);

    private:
        py::object mStructure;
};


std::vector<std::complex<double>> calculate_structure_factors_TAAM(py::object structure, const double d);

std::vector<std::complex<double>> calculate_structure_factors_IAM(py::object structure, const double d);
