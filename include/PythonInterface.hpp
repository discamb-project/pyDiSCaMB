#pragma once

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <complex>
#include <vector>

#include "DiscambStructureFactorCalculator.hpp"
#include "discamb/Scattering/SfCalculator.h"
#include "json.hpp"
#include "pybind11/complex.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11_json/pybind11_json.hpp"

namespace py = pybind11;

class PythonInterface : public DiscambStructureFactorCalculator {
   public:
    PythonInterface(py::object &structure, py::dict kwargs);

    virtual void set_indices(py::object &indices);
    virtual void set_d_min(const double d_min);
    virtual void update_structure(py::object &structure);

   private:
    py::object mStructure;
};
