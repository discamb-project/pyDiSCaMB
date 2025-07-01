#pragma once

#include <complex>
#include <string>
#include <vector>

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/StructuralParametersConverter.h"
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/Scattering/SfCalculator.h"
#include "json.hpp"

struct FCalcDerivatives;

class DiscambStructureFactorCalculator {
   public:
    DiscambStructureFactorCalculator(
        nlohmann::json calculator_parameters, discamb::Crystal crystal,
        std::vector<std::complex<double>> anomalous);
    virtual void init();
    ~DiscambStructureFactorCalculator();

    virtual std::vector<std::complex<double>> f_calc();

    virtual std::vector<FCalcDerivatives> d_f_calc_d_params();
    FCalcDerivatives d_f_calc_hkl_d_params(int h, int k, int l);

    virtual std::vector<discamb::TargetFunctionAtomicParamDerivatives>
    d_target_d_params(std::vector<std::complex<double>> d_target_d_f_calc);

    virtual std::vector<discamb::TargetFunctionAtomicParamDerivatives>
    selected_d_target_d_params(
        std::vector<std::complex<double>> d_target_d_f_calc, bool site,
        bool adp, bool occupancy, bool fp);

    std::vector<discamb::Vector3i> hkl;

    discamb::Crystal crystal;
    std::vector<std::complex<double>> anomalous;
    virtual void update_calculator();
    virtual discamb::SfCalculator *init_calc(
        discamb::Crystal crystal, nlohmann::json calculator_parameters);

   private:
    nlohmann::json mParams;
    discamb::SfCalculator *mCalculator;
    discamb::StructuralParametersConverter mConverter;
    std::vector<std::complex<double>> mFcalc;
    bool mStale;
};

struct FCalcDerivatives : discamb::SfDerivativesAtHkl {
    std::vector<int> hkl;
    std::complex<double> structure_factor;
    std::complex<double> fpDerivative{0.0, 0.0};
    std::complex<double> fdpDerivative{0.0, 0.0};

    // Translate discamb::Vector3 to std::vector
    std::vector<std::vector<std::complex<double>>> siteDerivatives() const;
};
