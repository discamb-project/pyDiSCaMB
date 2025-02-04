#include "DiscambStructureFactorCalculator.hpp"
#include "atom_assignment.hpp"
#include "crystal_conversion.hpp"

#include "json.hpp"

using namespace discamb;
using namespace std;


vector<vector<complex<double>>> FCalcDerivatives::siteDerivatives() const{
    vector<vector<complex<double>>> out;
    for (int i = 0; i < atomicPostionDerivatives.size(); i++){
        out.push_back({
            atomicPostionDerivatives[i].x,
            atomicPostionDerivatives[i].y,
            atomicPostionDerivatives[i].z
        });
    }
    return out;
}

vector<complex<double>> DiscambStructureFactorCalculator::f_calc(){
    // py::print("Discamb f_calc");
    vector<complex<double>> sf;
    sf.resize(hkl.size());
    vector<bool> count_atom_contribution (mCrystal.atoms.size(), true);
    mCalculator->calculateStructureFactors(mCrystal.atoms, hkl, sf, count_atom_contribution);
    return sf;
}

vector<FCalcDerivatives> DiscambStructureFactorCalculator::d_f_calc_d_params(){
    vector<FCalcDerivatives> out;
    out.resize(hkl.size());

    for (int i = 0; i < hkl.size(); i++){
        out[i] = d_f_calc_hkl_d_params(hkl[i].x, hkl[i].y, hkl[i].z);
    }
    return out;
}

FCalcDerivatives DiscambStructureFactorCalculator::d_f_calc_hkl_d_params(int h, int k, int l){
    FCalcDerivatives out;
    out.hkl = {h, k, l};
    vector<bool> count_atom_contribution (mCrystal.atoms.size(), true);
    mCalculator->calculateStructureFactorsAndDerivatives(
            out.hkl,
            out.structure_factor,
            out,
            count_atom_contribution
        );
    return out;
}

vector<TargetFunctionAtomicParamDerivatives> DiscambStructureFactorCalculator::d_target_d_params(vector<complex<double>> d_target_d_f_calc){
    assert(hkl.size() == d_target_d_f_calc.size());
    vector<complex<double>> sf;
    vector<TargetFunctionAtomicParamDerivatives> out;
    out.resize(mCrystal.atoms.size());
    vector<bool> count_atom_contribution( mCrystal.atoms.size(), true );

    /* TODO

    // Ensure correct convention
    structural_parameters_convention::AdpConvention ac;
    structural_parameters_convention::XyzCoordinateSystem xyzc;
    mCalculator->getParametersConvention(xyzc, ac);

    mCalculator->setParametersConvention(
        structural_parameters_convention::XyzCoordinateSystem::cartesian,
        structural_parameters_convention::AdpConvention::U_cart
    ),

    mCalculator->calculateStructureFactorsAndDerivatives(
        hkl,
        sf,
        out,
        d_target_d_f_calc,
        count_atom_contribution
    );
    mCalculator->setParametersConvention(xyzc, ac);
    return out;
    */
    mCalculator->calculateStructureFactorsAndDerivatives(
        mCrystal.atoms,
        hkl,
        sf,
        out,
        d_target_d_f_calc,
        count_atom_contribution
    );
    return out;
}
