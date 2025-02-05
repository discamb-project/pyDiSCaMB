#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/Scattering/SfCalculator.h"

#include <string>
#include <vector>
#include <complex>


struct FCalcDerivatives : discamb::SfDerivativesAtHkl {
    std::vector<int> hkl;
    std::complex<double> structure_factor;
    std::complex<double> fpDerivative {0.0, 0.0};
    std::complex<double> fdpDerivative {0.0, 0.0};

    // Translate discamb::Vector3 to std::vector
    std::vector<std::vector<std::complex<double>>> siteDerivatives() const;
};

class DiscambStructureFactorCalculator {
    public:
        DiscambStructureFactorCalculator() = default;
        DiscambStructureFactorCalculator(
            discamb::SfCalculator *calculator, 
            discamb::Crystal crystal, 
            std::vector<std::complex<double>> anomalous
        );
        // ~DiscambStructureFactorCalculator(); // TODO

        std::vector<std::complex<double>> f_calc();

        std::vector<FCalcDerivatives> d_f_calc_d_params();
        FCalcDerivatives d_f_calc_hkl_d_params(int h, int k, int l);
        std::vector<discamb::TargetFunctionAtomicParamDerivatives> d_target_d_params(std::vector<std::complex<double>> d_target_d_f_calc);
        
        std::vector<discamb::Vector3i> hkl;

    private:
        discamb::SfCalculator *mCalculator; // Pointer since abstract class
        discamb::Crystal mCrystal;
        std::vector<std::complex<double>> mAnomalous;
        void update_calculator();
};
