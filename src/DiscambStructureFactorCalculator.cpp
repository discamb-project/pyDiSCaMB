#include "DiscambStructureFactorCalculator.hpp"
#include "atom_assignment.hpp"

#include "discamb/CrystalStructure/StructuralParametersConverter.h"

#include "assert.hpp"

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


DiscambStructureFactorCalculator::DiscambStructureFactorCalculator(
    nlohmann::json calculator_parameters, 
    Crystal crystal, 
    vector<complex<double>> anomalous
) : 
    mCalculator(SfCalculator::create(crystal, calculator_parameters)), 
    mCrystal(crystal), 
    mAnomalous(anomalous),
    mConverter(mCrystal.unitCell)
{
    assert(mCrystal.atoms.size() > 0);
    assert(mAnomalous.size() > 0);
    assert(mCrystal.atoms.size() == mAnomalous.size());
    update_calculator();
}

vector<complex<double>> DiscambStructureFactorCalculator::f_calc(){
    update_calculator();
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
    update_calculator();
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
    update_calculator();
    assert(hkl.size() == d_target_d_f_calc.size());
    vector<complex<double>> sf;
    vector<TargetFunctionAtomicParamDerivatives> out;
    out.resize(mCrystal.atoms.size());
    vector<bool> count_atom_contribution( mCrystal.atoms.size(), true );

    mCalculator->calculateStructureFactorsAndDerivatives(
        mCrystal.atoms,
        hkl,
        sf,
        out,
        d_target_d_f_calc,
        count_atom_contribution
    );

    // Ensure correct convention (U_cart and Cartesian)
    structural_parameters_convention::AdpConvention ac = mCrystal.adpConvention;
    structural_parameters_convention::XyzCoordinateSystem xyzc = mCrystal.xyzCoordinateSystem;

    int idx, nAtoms = out.size();
    vector<complex<double> > adpIn(6), adpOut(6);
    Vector3<complex<double> > xyzIn, xyzOut;
    int i;
    // ADPs
    for (idx = 0; idx < nAtoms; idx++){
        if (out[idx].adp_derivatives.size() == 6){
            for (i = 0; i < 6; i++)
                adpIn[i] = out[idx].adp_derivatives[i];
            mConverter.convertDerivativesADP(adpIn, adpOut, ac, structural_parameters_convention::AdpConvention::U_cart);
            for (i = 0; i<6; i++)
                out[idx].adp_derivatives[i] = adpOut[i].real();
        }
    }
    // Coordinates
    for (idx = 0; idx<nAtoms; idx++){
        for (i = 0; i < 3; i++)
            xyzIn[i] = out[idx].atomic_position_derivatives[i];
        mConverter.convertDerivativesXyz(xyzIn, xyzOut, xyzc,  structural_parameters_convention::XyzCoordinateSystem::cartesian);
        for (i = 0; i < 3; i++)
            out[idx].atomic_position_derivatives[i] = xyzOut[i].real();
    }
    
    return out;
}

void DiscambStructureFactorCalculator::update_calculator(){
    // mCalculator->update(mCrystal.atoms); // Already handled since we pass atoms to calculations
    assert(mAnomalous.size() == mCrystal.atoms.size());
    mCalculator->setAnomalous(mAnomalous);
}
