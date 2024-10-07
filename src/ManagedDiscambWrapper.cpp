#include "ManagedDiscambWrapper.hpp"

#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/Scattering/taam_utilities.h"

#include "discamb/IO/MATTS_BankReader.h"

#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"

#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"

namespace py = pybind11;

using namespace std;
using namespace discamb;

vector<complex<double>> ManagedDiscambWrapper::f_calc(){
    update_atoms(mManager.mCrystal);
    return mManager.f_calc();
}


ManagedDiscambWrapper::FCalcManager ManagedDiscambWrapper::manager_setup(double d_min, FCalcMethod method){
    Crystal crystal;
    get_crystal(crystal);
    vector<Vector3i> hkl;
    get_hkl(d_min, hkl);
    string table = get_discamb_table_string();
    switch (method)
    {
    case FCalcMethod::IAM:{
        AnyScattererStructureFactorCalculator calculator = get_IAM_calculator(crystal, table);
        return FCalcManager(calculator, crystal, hkl);
        break;
    }
    case FCalcMethod::TAAM:{
        AnyScattererStructureFactorCalculator calculator = get_TAAM_calculator(crystal);
        return FCalcManager(calculator, crystal, hkl);
        break;
    }
    default:
        on_error::throwException("Invalid method requested");
        break;
    }
}

vector<complex<double>> ManagedDiscambWrapper::FCalcManager::f_calc(){
    vector<complex<double>> sf;
    mCalculator.update(mCrystal.atoms);
    mCalculator.calculateStructureFactors(mHkl, sf);
    return sf;
}
