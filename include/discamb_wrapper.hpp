#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/MathUtilities/Vector3.h"


#include <string>
#include <vector>
#include <complex>


std::string get_discamb_version();

void calculateSfTaamMinimal(
    const discamb::Crystal& crystal,
    const std::vector<discamb::Vector3i>& hkl,
    std::vector< std::complex<double> > &structureFactors);

// void cctbx_structure_to_discamb_crystal(const py::object structure, discamb::Crystal &crystal);
// std::string print_pdb_from_structure(const py::object structure);
