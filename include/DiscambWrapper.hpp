#pragma once

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/MathUtilities/Vector3.h"


#include <string>
#include <vector>
#include <complex>


// TODO: Declare cctbx to discamb funcs for crystal
// TODO: write throughput tests
// TODO: avoid recalculating hkl every time

namespace py = pybind11;


enum FCalcMethod {
    IAM,
    TAAM
};


class DiscambWrapper {
    public:
        DiscambWrapper(py::object structure, std::string table) : mStructure(structure), mCrystal(), mTableString(table) { structure_to_crystal(mStructure, mCrystal); };
        DiscambWrapper(py::object structure) { DiscambWrapper(structure, "electron-IT"); }
        
        void update_structure(const py::object structure);
        void update_crystal(const discamb::Crystal &crystal);

        std::vector<std::complex<double>> f_calc(const double d_min, FCalcMethod method);

    private:
        py::object mStructure;
        discamb::Crystal mCrystal;
        std::string mTableString;

        std::vector<std::complex<double>> f_calc_hkl(const std::vector<discamb::Vector3i> &hkl, FCalcMethod method);
        std::vector<discamb::Vector3i> get_hkl(const double d_min) const;

        static void structure_to_crystal(const py::object structure, discamb::Crystal &crystal);
        static void structure_to_hkl(const py::object structure, double d, std::vector<discamb::Vector3i> &hkl);
        bool is_crystal_valid(const discamb::Crystal &crystal) const;

};

// TEMPORARY
void calculateSfTaamMinimal(
    const discamb::Crystal& crystal,
    const std::vector<discamb::Vector3i>& hkl,
    std::vector< std::complex<double> > &structureFactors);

void calculateSfIamMinimal(
    const discamb::Crystal& crystal,
    const std::vector<discamb::Vector3i>& hkl,
    std::vector< std::complex<double> >& structureFactors);
