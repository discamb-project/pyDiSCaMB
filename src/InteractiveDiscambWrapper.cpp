#include "InteractiveDiscambWrapper.hpp"

#include "discamb/Scattering/IamFormFactorCalculationsManager.h"

namespace py = pybind11;

using namespace std;
using namespace discamb;

vector<complex<double>> InteractiveDiscambWrapper::f_calc(){
    update_atoms(mManager.mCrystal);
    return mManager.f_calc();
}

InteractiveDiscambWrapper::FCalcManager InteractiveDiscambWrapper::manager_setup(double d_min, string method){
    Crystal crystal;
    get_crystal(crystal);
    vector<Vector3i> hkl;
    get_hkl(d_min, hkl);
    shared_ptr<AtomicFormFactorCalculationsManager> formfactorCalculator;
    formfactorCalculator = shared_ptr<AtomicFormFactorCalculationsManager>(
        new IamFormFactorCalculationsManager(crystal, mTableString));
    AnyScattererStructureFactorCalculator calculator(crystal);
    calculator.setAtomicFormfactorManager(formfactorCalculator);
    return FCalcManager(calculator, crystal, hkl);
}

vector<complex<double>> InteractiveDiscambWrapper::FCalcManager::f_calc(){
    vector<complex<double>> sf;
    mCalculator.update(mCrystal.atoms);
    mCalculator.calculateStructureFactors(mHkl, sf);
    return sf;
}
