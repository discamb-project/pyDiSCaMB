#pragma once

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/Scattering/SfCalculator.h"

#include <string>
#include <vector>
#include <complex>

#include "DiscambStructureFactorCalculator.hpp"

namespace py = pybind11;

enum FCalcMethod {
    IAM,
    TAAM
};

class DiscambWrapper {
    public:
        DiscambWrapper(py::object structure, FCalcMethod method = FCalcMethod::IAM);

        DiscambWrapper from_TAAM_parameters(
            py::object structure,
            bool convert_to_electron_scattering,
            std::string bank_filepath,
            std::string assignment_log_filepath,
            std::string parameter_log_filepath,
            std::string multipolar_cif_output_filepath,
            double unit_cell_charge,
            bool perform_parameter_scaling_from_unit_cell_charge
        );

        void set_indices(py::object indices);
        void set_d_min(const double d_min);

        std::vector<std::complex<double>> f_calc();
        std::vector<std::complex<double>> f_calc(const double d_min);

        std::vector<FCalcDerivatives> d_f_calc_d_params();
        FCalcDerivatives d_f_calc_hkl_d_params(py::tuple hkl);
        FCalcDerivatives d_f_calc_hkl_d_params(int h, int k, int l);
        std::vector<discamb::TargetFunctionAtomicParamDerivatives> d_target_d_params(std::vector<std::complex<double>> d_target_d_f_calc);
        
    private:
        py::object mStructure;
        DiscambStructureFactorCalculator mDiscambCalculator;
};

std::vector<std::complex<double>> calculate_structure_factors_TAAM(py::object structure, const double d);

std::vector<std::complex<double>> calculate_structure_factors_IAM(py::object structure, const double d);
