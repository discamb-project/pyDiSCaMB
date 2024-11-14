#pragma once

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"

#include <string>
#include <vector>
#include <complex>
#include <utility>


namespace py = pybind11;


enum FCalcMethod {
    IAM,
    TAAM
};

struct FCalcDerivatives : discamb::SfDerivativesAtHkl {
    std::vector<int> hkl;
    std::complex<double> structure_factor;
    std::complex<double> fpDerivative {0.0, 0.0};
    std::complex<double> fdpDerivative {0.0, 0.0};

    // Translate discamb::Vector3 to std::vector
    std::vector<std::vector<std::complex<double>>> siteDerivatives();
};


class DiscambWrapper {
    public:
        DiscambWrapper(py::object structure, FCalcMethod method = FCalcMethod::IAM);
        void set_indices(py::object indices);
        void set_d_min(const double d_min);
        std::vector<std::complex<double>> f_calc();
        std::vector<std::complex<double>> f_calc(const double d_min);
        std::vector<FCalcDerivatives> d_f_calc_d_params();
        
    protected:
        py::object mStructure;
        discamb::Crystal mCrystal;
        discamb::AnyScattererStructureFactorCalculator mCalculator;
        std::vector<std::complex<double>> mAnomalous;
        std::vector<discamb::Vector3i> mHkl;


        void update();
        void init_crystal();
        discamb::AnyScattererStructureFactorCalculator get_calculator();
        void update_atoms();
        std::string get_discamb_table_string();
};

void set_IAM_calculator( discamb::AnyScattererStructureFactorCalculator &calculator, discamb::Crystal &crystal, std::string &table);
void set_TAAM_calculator( discamb::AnyScattererStructureFactorCalculator &calculator, discamb::Crystal &crystal);
