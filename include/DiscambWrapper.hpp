#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/MathUtilities/Vector3.h"


#include <string>
#include <vector>
#include <complex>


// TODO: Declare cctbx to discamb funcs for crystal
// TODO: write throughput tests
// TODO: avoid recalculating hkl every time

enum FCalcMethod {
    IAM,
    TAAM
};

class DiscambWrapper {
    public:
        DiscambWrapper() = default;
        DiscambWrapper(discamb::Crystal &crystal, std::string table) : mCrystal(crystal), mTableString(table) {};
        DiscambWrapper(discamb::Crystal &crystal) { DiscambWrapper(crystal, "electron-IT"); }
        
        void update_crystal(const discamb::Crystal &crystal);

        std::vector<std::complex<double>> f_calc_hkl(const std::vector<discamb::Vector3i> &hkl, FCalcMethod method);

        std::vector<std::complex<double>> f_calc(const double d_min, FCalcMethod method);

    private:
        discamb::Crystal mCrystal;
        std::string mTableString;
        virtual std::vector<discamb::Vector3i> get_hkl(const double d_min) const {};
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
