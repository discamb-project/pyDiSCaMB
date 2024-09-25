#pragma once

#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"

#include "DiscambWrapper.hpp"

namespace py = pybind11;

class InteractiveDiscambWrapper : public DiscambWrapper {
    using DiscambWrapper::DiscambWrapper;

    public:
        InteractiveDiscambWrapper(
            py::object structure, 
            double d_min, 
            FCalcMethod method
        ): 
            DiscambWrapper(structure),
            mManager(manager_setup(d_min, method))
        {};

        std::vector<std::complex<double>> f_calc();

    private:
        
        // a little context manager to avoid some data transfer and re-calculations
        class FCalcManager {

            public:
                FCalcManager(
                    discamb::AnyScattererStructureFactorCalculator &calculator, 
                    discamb::Crystal &crystal,
                    std::vector<discamb::Vector3i> &hkl
                ):
                    mCalculator(calculator),
                    mCrystal(crystal),
                    mHkl(hkl)
                {};
                
                std::vector<std::complex<double>> f_calc();

            private:
                discamb::AnyScattererStructureFactorCalculator mCalculator;
                discamb::Crystal mCrystal;
                std::vector<discamb::Vector3i> mHkl;
            
            friend class InteractiveDiscambWrapper;

        };
        FCalcManager manager_setup(double d_min, FCalcMethod method);
        FCalcManager mManager;
        
};
