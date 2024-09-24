#include "DiscambWrapper.hpp"

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


// Forward declarations
void calculateSfTaamMinimal(
    const Crystal& crystal,
    const vector<Vector3i>& hkl,
    vector< complex<double> > &structureFactors);

void calculateSfIamMinimal(
    const Crystal& crystal,
    const vector<Vector3i>& hkl,
    vector< complex<double> >& structureFactors);

void DiscambWrapper::update_structure(const py::object structure){
    Crystal new_crystal;
    structure_to_crystal(structure, new_crystal);
    update_crystal(new_crystal);
    mStructure = structure;
}

void DiscambWrapper::update_crystal(const Crystal &crystal){
    if (!is_crystal_valid(crystal)){
        on_error::throwException("Incompatible crystal");
    }
    for (int i = 0; i < crystal.atoms.size(); i++){
        mCrystal.atoms[i].coordinates[0] = crystal.atoms[i].coordinates[0];
        mCrystal.atoms[i].coordinates[1] = crystal.atoms[i].coordinates[1];
        mCrystal.atoms[i].coordinates[2] = crystal.atoms[i].coordinates[2];
        for (int j = 0; j < crystal.atoms[i].adp.size(); j++){
            mCrystal.atoms[i].adp[j] = crystal.atoms[i].adp[j];
        }
        mCrystal.atoms[i].occupancy = crystal.atoms[i].occupancy;
    }
}
    

vector<complex<double>> DiscambWrapper::f_calc_hkl(const vector<Vector3i> &hkl, FCalcMethod method){
    vector<complex<double>> sf;
    sf.reserve(hkl.size());
    switch (method)
    {
        case FCalcMethod::IAM:
            calculateSfIamMinimal(mCrystal, hkl, sf);
            break;
        case FCalcMethod::TAAM:
            calculateSfTaamMinimal(mCrystal, hkl, sf);
            break;
        default:
            on_error::throwException("Invalid method requested");
            break;
    }
    return sf;
}

vector<complex<double>> DiscambWrapper::f_calc(const double d_min, FCalcMethod method){
    vector<Vector3i> hkl = get_hkl(d_min);
    return f_calc_hkl(hkl, method);
}

vector<Vector3i> DiscambWrapper::get_hkl(const double d_min) const{
    vector<Vector3i> hkl;
    structure_to_hkl(mStructure, d_min, hkl);
    return hkl;
}

bool DiscambWrapper::is_crystal_valid(const Crystal &crystal) const{
    // TODO
    return true;
}
void DiscambWrapper::structure_to_crystal(const py::object structure, Crystal &crystal){
    py::tuple cell_params = structure.attr("unit_cell")().attr("parameters")();
    UnitCell cell {
        cell_params[0].cast<double>(), 
        cell_params[1].cast<double>(), 
        cell_params[2].cast<double>(), 
        cell_params[3].cast<double>(), 
        cell_params[4].cast<double>(), 
        cell_params[5].cast<double>()
        };
    crystal.unitCell = cell;

    vector<SpaceGroupOperation> symops;
    for (auto symop_py : structure.attr("crystal_symmetry")().attr("space_group")()){
        string symop_str = symop_py.attr("as_xyz")().cast<string>();
        SpaceGroupOperation symop {symop_str};
        symops.push_back(symop);
    }
    SpaceGroup space_group;
    space_group.set(symops);
    crystal.spaceGroup = space_group;

    vector<AtomInCrystal> atoms;
    for (auto scatterer_py : structure.attr("scatterers")()){
        AtomInCrystal atom;

        py::tuple xyz_py = scatterer_py.attr("site");
        atom.coordinates[0] = xyz_py[0].cast<float>();
        atom.coordinates[1] = xyz_py[1].cast<float>();
        atom.coordinates[2] = xyz_py[2].cast<float>();

        atom.coordinates_sigma[0] = 0.0;
        atom.coordinates_sigma[1] = 0.0;
        atom.coordinates_sigma[2] = 0.0;

        atom.adp.clear();
        atom.adp.push_back(scatterer_py.attr("u_iso").cast<float>()); // U11
        // atom.adp.push_back(scatterer_py.attr("u_iso").cast<float>()); // U22
        // atom.adp.push_back(scatterer_py.attr("u_iso").cast<float>()); // U33
        // atom.adp.push_back(0.0); // U12
        // atom.adp.push_back(0.0); // U13
        // atom.adp.push_back(0.0); // U23
        atom.adp_sigma.clear();
        atom.adp_sigma.push_back(0.0);
        // for(int i = 0; i < 6; i++){atom.adp_sigma.push_back(0.0);}
        
        atom.label = scatterer_py.attr("scattering_type").cast<string>();

        atom.occupancy = scatterer_py.attr("occupancy").cast<float>();
        atom.multiplicity = scatterer_py.attr("multiplicity")().cast<float>();
        atom.occupancy_sigma = 0.0;

        // This seems to be unused anyway
        atom.siteSymetry.clear();
        SpaceGroupOperation identity = SpaceGroupOperation();
        atom.siteSymetry.push_back(identity);

        atoms.push_back(atom);
    }
    crystal.atoms = atoms;

    crystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;
    crystal.adpConvention = structural_parameters_convention::AdpConvention::U_cif;
}

void DiscambWrapper::structure_to_hkl(const py::object structure, double d, vector<Vector3i> &hkl){
    // assume anomolous_flag is False
    py::object miller_py = structure.attr("build_miller_set")(false, d);

    hkl.clear();
    for (auto hkl_py_auto : miller_py.attr("indices")().attr("as_vec3_double")()){
        py::tuple hkl_py = hkl_py_auto.cast<py::tuple>();
        Vector3i hkl_i {
            static_cast<int>(hkl_py[0].cast<double>()),
            static_cast<int>(hkl_py[1].cast<double>()),
            static_cast<int>(hkl_py[2].cast<double>())
        };
        hkl.push_back(hkl_i);
    }
}


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
