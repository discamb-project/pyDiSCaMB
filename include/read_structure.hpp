#pragma once

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <complex>
#include <string>
#include <vector>

#include "discamb/CrystalStructure/Crystal.h"
#include "pybind11/complex.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

discamb::Crystal crystal_from_xray_structure(const py::object &structure);

std::vector<std::complex<double>> anomalous_from_xray_structure(
    const py::object &structure);

std::string table_from_xray_structure(const py::object &structure);

void update_crystal_from_xray_structure(discamb::Crystal &crystal,
                                        const py::object &structure);

void update_anomalous_from_xray_structure(
    std::vector<std::complex<double>> &anomalous, const py::object &structure);
