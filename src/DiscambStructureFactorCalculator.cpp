#include "DiscambStructureFactorCalculator.hpp"

#include <chrono>

#include "assert.hpp"
#include "atom_assignment.hpp"
#include "discamb/CrystalStructure/StructuralParametersConverter.h"

using namespace discamb;
using namespace std;

DiscambStructureFactorCalculator::DiscambStructureFactorCalculator(
    nlohmann::json calculator_parameters, Crystal crystal,
    vector<complex<double>> anomalous)
    : mParams(calculator_parameters),
      mCalculator(),
      crystal(crystal),
      anomalous(anomalous),
      mConverter(crystal.unitCell),
      mStale(true),
      mFcalc({{0.0, 0.0}}) {
    // init();
}

void DiscambStructureFactorCalculator::init() {
    auto runtime_start = chrono::high_resolution_clock::now();
    mCalculator = init_calc(crystal, mParams);
    auto runtime_end = chrono::high_resolution_clock::now();
    chrono::duration<double> delta =
        chrono::duration_cast<chrono::nanoseconds>(runtime_end - runtime_start);
    double time = delta.count();
    cout << "init_calc: " << time << endl;
    assert(crystal.atoms.size() > 0);
    assert(anomalous.size() > 0);
    assert(crystal.atoms.size() == anomalous.size());
    update_calculator();
    runtime_end = chrono::high_resolution_clock::now();
    delta =
        chrono::duration_cast<chrono::nanoseconds>(runtime_end - runtime_start);
    time = delta.count();
    cout << "init: " << time << endl;
}
DiscambStructureFactorCalculator::~DiscambStructureFactorCalculator() {
    delete mCalculator;
}

vector<complex<double>> DiscambStructureFactorCalculator::f_calc() {
    if (mStale) {
        mFcalc.resize(hkl.size());
        vector<bool> count_atom_contribution(crystal.atoms.size(), true);
        auto runtime_start = chrono::high_resolution_clock::now();
        mCalculator->calculateStructureFactors(
            crystal.atoms, hkl, mFcalc, count_atom_contribution);
        auto runtime_end = chrono::high_resolution_clock::now();
        chrono::duration<double> delta =
            chrono::duration_cast<chrono::nanoseconds>(runtime_end -
                                                       runtime_start);
        double time = delta.count();
        cout << "f_calc: " << time << endl;
        mStale = false;
    }
    return mFcalc;
}

vector<FCalcDerivatives> DiscambStructureFactorCalculator::d_f_calc_d_params() {
    vector<FCalcDerivatives> out;
    out.resize(hkl.size());

    for (int i = 0; i < hkl.size(); i++) {
        out[i] = d_f_calc_hkl_d_params(hkl[i].x, hkl[i].y, hkl[i].z);
    }
    return out;
}

FCalcDerivatives DiscambStructureFactorCalculator::d_f_calc_hkl_d_params(
    int h, int k, int l) {
    FCalcDerivatives out;
    out.hkl = {h, k, l};
    vector<bool> count_atom_contribution(crystal.atoms.size(), true);
    mCalculator->calculateStructureFactorsAndDerivatives(
        out.hkl, out.structure_factor, out, count_atom_contribution);
    return out;
}

vector<TargetFunctionAtomicParamDerivatives>
DiscambStructureFactorCalculator::d_target_d_params(
    vector<complex<double>> d_target_d_f_calc) {
    return selected_d_target_d_params(
        d_target_d_f_calc, true, true, true, true);
}

vector<TargetFunctionAtomicParamDerivatives>
DiscambStructureFactorCalculator::selected_d_target_d_params(
    vector<complex<double>> d_target_d_f_calc, bool site, bool adp,
    bool occupancy, bool fp) {
    assert(hkl.size() == d_target_d_f_calc.size());
    mFcalc.resize(hkl.size());
    vector<TargetFunctionAtomicParamDerivatives> out;
    out.resize(crystal.atoms.size());
    vector<bool> count_atom_contribution(crystal.atoms.size(), true);

    DerivativesSelector ds{site, adp, occupancy, fp};
    auto runtime_start = chrono::high_resolution_clock::now();
    mCalculator->calculateStructureFactorsAndDerivatives(
        crystal.atoms,
        hkl,
        mFcalc,
        out,
        d_target_d_f_calc,
        count_atom_contribution,
        ds);
    auto runtime_end = chrono::high_resolution_clock::now();
    chrono::duration<double> delta =
        chrono::duration_cast<chrono::nanoseconds>(runtime_end - runtime_start);
    double time = delta.count();
    cout << "grad: " << time << endl;

    mStale = false;

    // Ensure correct convention (U_cart and Cartesian)
    structural_parameters_convention::AdpConvention ac = crystal.adpConvention;
    structural_parameters_convention::XyzCoordinateSystem xyzc =
        crystal.xyzCoordinateSystem;

    int idx, nAtoms = out.size();
    vector<complex<double>> adpIn(6), adpOut(6);
    Vector3<complex<double>> xyzIn, xyzOut;
    int i;
    // ADPs
    for (idx = 0; idx < nAtoms; idx++) {
        if (out[idx].adp_derivatives.size() == 6) {
            for (i = 0; i < 6; i++) adpIn[i] = out[idx].adp_derivatives[i];
            mConverter.convertDerivativesADP(
                adpIn,
                adpOut,
                ac,
                structural_parameters_convention::AdpConvention::U_cart);
            for (i = 0; i < 6; i++)
                out[idx].adp_derivatives[i] = adpOut[i].real();
        }
    }
    // Coordinates
    for (idx = 0; idx < nAtoms; idx++) {
        for (i = 0; i < 3; i++)
            xyzIn[i] = out[idx].atomic_position_derivatives[i];
        mConverter.convertDerivativesXyz(
            xyzIn,
            xyzOut,
            xyzc,
            structural_parameters_convention::XyzCoordinateSystem::cartesian);
        for (i = 0; i < 3; i++)
            out[idx].atomic_position_derivatives[i] = xyzOut[i].real();
    }

    return out;
}

void DiscambStructureFactorCalculator::update_calculator() {
    // mCalculator->update(crystal.atoms); // Already handled since we pass
    // atoms to calculations
    assert(anomalous.size() == crystal.atoms.size());
    auto runtime_start = chrono::high_resolution_clock::now();
    mCalculator->setAnomalous(anomalous);
    mConverter.set(crystal.unitCell);
    auto runtime_end = chrono::high_resolution_clock::now();
    chrono::duration<double> delta =
        chrono::duration_cast<chrono::nanoseconds>(runtime_end - runtime_start);
    double time = delta.count();
    cout << "update calculator: " << time << endl;
    mStale = true;
}

discamb::SfCalculator* DiscambStructureFactorCalculator::init_calc(
    discamb::Crystal crystal, nlohmann::json calculator_parameters) {
    return SfCalculator::create(crystal, calculator_parameters);
}

vector<vector<complex<double>>> FCalcDerivatives::siteDerivatives() const {
    vector<vector<complex<double>>> out;
    for (int i = 0; i < atomicPostionDerivatives.size(); i++) {
        out.push_back({atomicPostionDerivatives[i].x,
                       atomicPostionDerivatives[i].y,
                       atomicPostionDerivatives[i].z});
    }
    return out;
}
