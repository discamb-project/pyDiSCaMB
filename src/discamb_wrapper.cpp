#include "discamb_wrapper.hpp"

#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"

#include "discamb/BasicUtilities/OnError.h"

#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"

#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/MATTS_BankReader.h"

#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/Scattering/MATTS_Default.h"
#include "discamb/Scattering/taam_utilities.h"


using namespace discamb;
using namespace std;


void calculateSfTaamMinimal(
    const Crystal& crystal,
    const vector<Vector3i>& hkl,
    vector< complex<double> > &structureFactors)
{

    bool electronScattering = true;
    bool scaleHcParameters = true;

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

    AnyScattererStructureFactorCalculator sfCalc(crystal);
    sfCalc.setAtomicFormfactorManager(formfactorCalculator);
    sfCalc.calculateStructureFactors(hkl, structureFactors);

}

void calculateSfIamMinimal(
    const Crystal& crystal,
    const vector<Vector3i>& hkl,
    vector< complex<double> >& structureFactors)
{
    bool electronScattering = true;
    //table can be "Waasmeier-Kirfel", "IT92", "electron-IT"
    string table = "electron-IT";
    AnyIamCalculator iamCalculator(crystal, electronScattering, table);
    vector<bool> countAtomContribution(crystal.atoms.size(), true);
    iamCalculator.calculateStructureFactors(crystal.atoms, hkl, structureFactors, countAtomContribution);
}

