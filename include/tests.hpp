#pragma once

#include <complex>
#include <string>
#include <vector>

#include "discamb/Scattering/SfCalculator.h"
#include "pybind11/pybind11.h"

namespace py = pybind11;

std::vector<std::complex<double>> f_calc_custom_gaussian_parameters(
    py::object structure, std::vector<std::vector<int>> hkl,
    std::vector<std::string> atom_labels, std::vector<std::vector<double>> a,
    std::vector<std::vector<double>> b, std::vector<double> c);
