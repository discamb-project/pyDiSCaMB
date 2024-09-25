#include "InteractiveDiscambWrapper.hpp"

namespace py = pybind11;

using namespace std;
using namespace discamb;

vector<complex<double>> InteractiveDiscambWrapper::f_calc(){
    update_atoms(mManager.mCrystal);
    return mManager.f_calc();
}

InteractiveDiscambWrapper::FCalcManager InteractiveDiscambWrapper::manager_setup(double d_min, string method){
    bool electron_scattering = mTableString == "electron-IT";
    Crystal crystal;
    get_crystal(crystal);
    vector<Vector3i> hkl;
    get_hkl(d_min, hkl);
    AnyIamCalculator calculator {crystal, electron_scattering, mTableString};
    return FCalcManager(calculator, crystal, hkl);
}

vector<complex<double>> InteractiveDiscambWrapper::FCalcManager::f_calc(){
    vector<complex<double>> sf;
    mCalculator.update(mCrystal.atoms);
    vector<bool> countAtomContribution(mCrystal.atoms.size(), true);
    mCalculator.calculateStructureFactors(mCrystal.atoms, mHkl, sf, countAtomContribution);
    return sf;
}
