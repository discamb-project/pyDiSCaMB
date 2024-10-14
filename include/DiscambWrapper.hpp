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


class DiscambWrapper {
    public:
        DiscambWrapper(py::object structure, FCalcMethod method = FCalcMethod::IAM);
        
        std::vector<std::complex<double>> f_calc(const double d_min);

        
    protected:
        py::object mStructure;
        discamb::Crystal mCrystal;
        discamb::AnyScattererStructureFactorCalculator mCalculator;
        std::vector<std::complex<double>> mAnomalous;


        void update();
        void init_crystal();
        void get_hkl(double d, std::vector<discamb::Vector3i> &hkl);
        discamb::AnyScattererStructureFactorCalculator get_calculator();
        void update_atoms();
        std::string get_discamb_table_string();
};

void set_IAM_calculator( discamb::AnyScattererStructureFactorCalculator &calculator, discamb::Crystal &crystal, std::string &table);
void set_TAAM_calculator( discamb::AnyScattererStructureFactorCalculator &calculator, discamb::Crystal &crystal);
