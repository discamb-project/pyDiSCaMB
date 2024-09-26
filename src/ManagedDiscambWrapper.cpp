#include "ManagedDiscambWrapper.hpp"

#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/Scattering/MATTS_Default.h"
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

// Forward declarations
AnyScattererStructureFactorCalculator get_IAM_calculator(Crystal &crystal, string &table);
AnyScattererStructureFactorCalculator get_TAAM_calculator(Crystal &crystal);


ManagedDiscambWrapper::FCalcManager ManagedDiscambWrapper::manager_setup(double d_min, FCalcMethod method){
    Crystal crystal;
    get_crystal(crystal);
    vector<Vector3i> hkl;
    get_hkl(d_min, hkl);
    switch (method)
    {
    case FCalcMethod::IAM:{
        AnyScattererStructureFactorCalculator calculator = get_IAM_calculator(crystal, mTableString);
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

AnyScattererStructureFactorCalculator get_IAM_calculator(Crystal &crystal, string &table){
    shared_ptr<AtomicFormFactorCalculationsManager> formfactorCalculator;
    formfactorCalculator = shared_ptr<AtomicFormFactorCalculationsManager>(
        new IamFormFactorCalculationsManager(crystal, table));
    AnyScattererStructureFactorCalculator calculator(crystal);
    calculator.setAtomicFormfactorManager(formfactorCalculator);
    return calculator;
}

AnyScattererStructureFactorCalculator get_TAAM_calculator(Crystal &crystal){

    // set bank parameters
    MATTS_BankReader bankReader;
    vector<AtomType> atomTypes;
    vector<AtomTypeHC_Parameters> hcParameters;
    BankSettings bankSettings;

    string bankString;
    stringstream bankStream;
    default_ubdb_bank_string(bankString);
    bankStream << bankString;
    bankReader.read(bankStream, atomTypes, hcParameters, bankSettings, true);

    // assign atom types

    CrystalAtomTypeAssigner assigner;
    assigner.setAtomTypes(atomTypes);
    assigner.setDescriptorsSettings(DescriptorsSettings());
    vector < LocalCoordinateSystem<AtomInCrystalID> > lcs;
    vector<int> types;
    assigner.assign(crystal, types, lcs);

    // get TAAM parameters with unit cell charge scaled/shifted to 0

    HC_ModelParameters multipoleModelPalameters;

    vector<int> atomicNumbers;
    vector<double> nuclearCharge;
    vector<int> nonMultipolarAtoms;
    crystal_structure_utilities::atomicNumbers(crystal, atomicNumbers);
    for (int z : atomicNumbers)
        nuclearCharge.push_back(z);

    vector<double> multiplicityTimesOccupancy;
    for (auto& atom : crystal.atoms)
        multiplicityTimesOccupancy.push_back(atom.multiplicity * atom.occupancy);
    double unitCellCharge = 0.0;
    taam_utilities::type_assignment_to_HC_parameters(
            hcParameters, types, multiplicityTimesOccupancy, atomicNumbers, unitCellCharge,
            multipoleModelPalameters, true, nonMultipolarAtoms);

    vector<shared_ptr<LocalCoordinateSystemInCrystal> > lcaCalculators;
    for (auto coordinateSystem : lcs)
        lcaCalculators.push_back(
            shared_ptr<LocalCoordinateSystemInCrystal>(
                new LocalCoordinateSystemCalculator(coordinateSystem, crystal)));


    //calculate electron structure factors 

    shared_ptr<AtomicFormFactorCalculationsManager> formfactorCalculator;

    shared_ptr<AtomicFormFactorCalculationsManager> hcManager =
        shared_ptr<AtomicFormFactorCalculationsManager>(
            new HcFormFactorCalculationsManager(crystal, multipoleModelPalameters, lcaCalculators));
        
    formfactorCalculator = shared_ptr<AtomicFormFactorCalculationsManager>(
        new ElectronFromXrayFormFactorCalculationManager(
            crystal.unitCell,
            nuclearCharge,
            hcManager));

    AnyScattererStructureFactorCalculator calculator(crystal);
    calculator.setAtomicFormfactorManager(formfactorCalculator);
    return calculator;
}
