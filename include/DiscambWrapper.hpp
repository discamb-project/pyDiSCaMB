#pragma once

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/MathUtilities/Vector3.h"


#include <string>
#include <vector>
#include <complex>
#include <utility>


// TODO: write throughput tests
// TODO: avoid recalculating hkl every time

namespace py = pybind11;


enum FCalcMethod {
    IAM,
    TAAM
};


class DiscambWrapper {
    public:
        DiscambWrapper(py::object structure, std::string table = "electron-IT") : 
            mStructure(std::move(structure)),
            mTableString(table) {};
        
        std::vector<std::complex<double>> f_calc(const double d_min, FCalcMethod method);
        std::vector<std::complex<double>> f_calc_IAM(const double d_min) { return f_calc(d_min, FCalcMethod::IAM); };
        std::vector<std::complex<double>> f_calc_TAAM(const double d_min) { return f_calc(d_min, FCalcMethod::TAAM); };

        void test_get_crystal();
    private:
        py::object mStructure;
        std::string mTableString;

        void f_calc_hkl(const std::vector<discamb::Vector3i> &hkl, FCalcMethod method, std::vector<std::complex<double>> &sf);

        void get_crystal(discamb::Crystal &crystal);
        void get_hkl(double d, std::vector<discamb::Vector3i> &hkl);

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
