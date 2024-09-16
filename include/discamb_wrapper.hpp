#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/MathUtilities/Vector3.h"


#include <string>
#include <vector>
#include <complex>


void calculateSfTaamMinimal(
    const discamb::Crystal& crystal,
    const std::vector<discamb::Vector3i>& hkl,
    std::vector< std::complex<double> > &structureFactors);

void calculateSfIamMinimal(
    const discamb::Crystal& crystal,
    const std::vector<discamb::Vector3i>& hkl,
    std::vector< std::complex<double> >& structureFactors);
