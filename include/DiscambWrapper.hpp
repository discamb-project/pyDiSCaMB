#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/MathUtilities/Vector3.h"


#include <string>
#include <vector>
#include <complex>


// TODO: Declare cctbx to discamb funcs for crystal
// TODO: write throughput tests

enum FCalcMethod {
    IAM,
    TAAM
};

class DiscambWrapper {
    DiscambWrapper(discamb::Crystal crystal, std::string table) : mCrystal(crystal), mTableString(table) [];
    DiscambWrapper(discamb::Crystal crystal) { DiscambWrapper(crystal, "electron-IT"); }
    
    // for this one: just do name lookup (ignore atom type, just ensure label is same)
    // also call is_crystal_valid
    void update_crystal(const discamb::Crystal crystal);

    std::vector<std::complex<double>> f_calc(const double d_min, FCalcMethod method);
    std::vector<std::complex<double>> f_calc(const double d_min) { f_calc(d_min, FCalcMethod::IAM); };

    private:
        discamb::Crystal mCrystal;
        std::string mTableString;
        std::vector<discamb::Vector3i> hkl(const double d_min) const;
        bool is_crystal_valid(const discamb::Crystal crystal) const;
};
