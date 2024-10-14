#include "DiscambWrapper.hpp"

#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"

#include "discamb/BasicUtilities/on_error.h"

#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"

#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/MATTS_BankReader.h"

#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/Scattering/taam_utilities.h"


using namespace discamb;
using namespace std;


DiscambWrapper::DiscambWrapper(py::object structure, FCalcMethod method) :
            mStructure(std::move(structure)),
            mAnomalous(std::vector<std::complex<double>> {}),
            mCrystal(),
            mCalculator(mCrystal)
            {
                init_crystal();
                mCalculator = discamb::AnyScattererStructureFactorCalculator(mCrystal);
                switch (method)
                {
                case FCalcMethod::IAM: {
                    string table = get_discamb_table_string();
                    set_IAM_calculator(mCalculator, mCrystal, table);
                    break;
                }
                case FCalcMethod::TAAM: {
                    set_TAAM_calculator(mCalculator, mCrystal);
                    break;
                }
                default:
                    break;
                }
            };

void DiscambWrapper::update(){
    update_atoms();
    mCalculator.setAnoumalous(mAnomalous);
    mCalculator.update(mCrystal.atoms);
}


vector<complex<double>> DiscambWrapper::f_calc(const double d_min){
    vector<Vector3i> hkl;
    get_hkl(d_min, hkl);
    vector<complex<double>> sf;
    sf.resize(hkl.size());
    update();
    mCalculator.calculateStructureFactors(hkl, sf);
    return sf;
}


void DiscambWrapper::init_crystal(){
    py::tuple cell_params = mStructure.attr("unit_cell")().attr("parameters")();
    mCrystal.unitCell.set(
        cell_params[0].cast<double>(), 
        cell_params[1].cast<double>(), 
        cell_params[2].cast<double>(), 
        cell_params[3].cast<double>(), 
        cell_params[4].cast<double>(), 
        cell_params[5].cast<double>()
    );

    vector<SpaceGroupOperation> symops;
    for (auto symop_py : mStructure.attr("crystal_symmetry")().attr("space_group")()){
        string symop_str = symop_py.attr("as_xyz")().cast<string>();
        SpaceGroupOperation symop {symop_str};
        symops.push_back(symop);
    }
    mCrystal.spaceGroup.set(symops);


    int num_atoms = mStructure.attr("scatterers")().attr("size")().cast<int>();
    mCrystal.atoms.clear();
    mCrystal.atoms.assign(num_atoms, AtomInCrystal());
    mAnomalous.clear();
    mAnomalous.assign(num_atoms, complex<double> (0.0, 0.0));
    update_atoms();
    

    mCrystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;
    // Assume all adps use the same convention
    if (mStructure.attr("use_u_iso")().attr("__getitem__")(0).cast<bool>()){
        mCrystal.adpConvention = structural_parameters_convention::AdpConvention::U_cif;
    }
    else{
        mCrystal.adpConvention = structural_parameters_convention::AdpConvention::U_star;
    }
}

void DiscambWrapper::get_hkl(double d, vector<Vector3i> &hkl){
    bool anomalous_flag = mStructure.attr("scatterers")().attr("count_anomalous")().cast<int>() != 0;
    py::object miller_py = mStructure.attr("build_miller_set")(anomalous_flag, d);

    for (auto hkl_py_auto : miller_py.attr("indices")().attr("as_vec3_double")()){
        py::tuple hkl_py = hkl_py_auto.cast<py::tuple>();
        hkl.push_back(Vector3i {
            static_cast<int>(hkl_py[0].cast<double>()),
            static_cast<int>(hkl_py[1].cast<double>()),
            static_cast<int>(hkl_py[2].cast<double>())
        });
    }
}

void DiscambWrapper::update_atoms(){
    int idx = 0;
    for (auto scatterer_py : mStructure.attr("scatterers")()){
        py::tuple xyz_py = scatterer_py.attr("site");
        mCrystal.atoms[idx].coordinates[0] = xyz_py[0].cast<float>();
        mCrystal.atoms[idx].coordinates[1] = xyz_py[1].cast<float>();
        mCrystal.atoms[idx].coordinates[2] = xyz_py[2].cast<float>();

        mCrystal.atoms[idx].coordinates_sigma[0] = 0.0;
        mCrystal.atoms[idx].coordinates_sigma[1] = 0.0;
        mCrystal.atoms[idx].coordinates_sigma[2] = 0.0;

        mCrystal.atoms[idx].coordinates_precision[0] = 0.0;
        mCrystal.atoms[idx].coordinates_precision[1] = 0.0;
        mCrystal.atoms[idx].coordinates_precision[2] = 0.0;

        mCrystal.atoms[idx].adp.clear();
        if (scatterer_py.attr("flags").attr("use_u_iso")().cast<bool>()){
            mCrystal.atoms[idx].adp.push_back(scatterer_py.attr("u_iso").cast<float>()); // U11
        }
        else {
            py::tuple u_star = scatterer_py.attr("u_star");
            mCrystal.atoms[idx].adp.push_back(u_star[0].cast<float>()); // U11
            mCrystal.atoms[idx].adp.push_back(u_star[1].cast<float>()); // U22
            mCrystal.atoms[idx].adp.push_back(u_star[2].cast<float>()); // U33
            mCrystal.atoms[idx].adp.push_back(u_star[3].cast<float>()); // U12
            mCrystal.atoms[idx].adp.push_back(u_star[4].cast<float>()); // U13
            mCrystal.atoms[idx].adp.push_back(u_star[5].cast<float>()); // U23
        }
        mCrystal.atoms[idx].adp_sigma.clear();
        mCrystal.atoms[idx].adp_precision.clear();
        for(int i = 0; i < mCrystal.atoms[idx].adp.size(); i++){
            mCrystal.atoms[idx].adp_sigma.push_back(0.0);
            mCrystal.atoms[idx].adp_precision.push_back(0.0);
        }
        
        mCrystal.atoms[idx].label = scatterer_py.attr("scattering_type").cast<string>();

        mCrystal.atoms[idx].occupancy = scatterer_py.attr("occupancy").cast<float>();
        mCrystal.atoms[idx].multiplicity = scatterer_py.attr("multiplicity")().cast<float>();
        mCrystal.atoms[idx].occupancy_sigma = 0.0;

        // This seems to be unused anyway
        mCrystal.atoms[idx].siteSymetry.clear();
        mCrystal.atoms[idx].siteSymetry.push_back(SpaceGroupOperation());

        mAnomalous[idx] = complex<double> {
            scatterer_py.attr("fp").cast<double>(),
            scatterer_py.attr("fdp").cast<double>()
        };
        idx++;
    }
}

string DiscambWrapper::get_discamb_table_string(){
    auto py_table_string = mStructure.attr("scattering_type_registry_params").attr("table");
    // default (i.e. no table provided)
    if (py::isinstance<py::none>(py_table_string)){
        return "Waasmeier-Kirfel";
    }
    string cctbx_table_string = py_table_string.cast<string>();
    string discamb_table_string;
    // Valid cctbx entries:
    // ["n_gaussian", "it1992", "wk1995", "xray", "electron", "neutron"]
    if (cctbx_table_string == "electron"){
        discamb_table_string = "electron-cctbx";
    }
    else if (cctbx_table_string == "it1992"){
        discamb_table_string = "IT92";
    }
    else if (cctbx_table_string == "wk1995"){
        discamb_table_string = "Waasmeier-Kirfel";
    }
    else {
        throw std::invalid_argument("Scattering table \"" + cctbx_table_string + "\" is not recognized.");
    }
    return discamb_table_string;
}

void set_IAM_calculator(AnyScattererStructureFactorCalculator &calculator, Crystal &crystal, string &table){
    shared_ptr<AtomicFormFactorCalculationsManager> formfactorCalculator;
    formfactorCalculator = shared_ptr<AtomicFormFactorCalculationsManager>(
        new IamFormFactorCalculationsManager(crystal, table));
    calculator.setAtomicFormfactorManager(formfactorCalculator);
}

void set_TAAM_calculator(AnyScattererStructureFactorCalculator &calculator, Crystal &crystal){

    // set bank parameters
    MATTS_BankReader bankReader;
    vector<AtomType> atomTypes;
    vector<AtomTypeHC_Parameters> hcParameters;
    BankSettings bankSettings;

    string TAAM_root = py::module::import("pydiscamb").attr("_get_TAAM_databank_directory")().cast<string>();

    // Try to load the MATTS bank
    ifstream bankStream { TAAM_root + "/MATTS2021databank.txt" };
    if (!bankStream.good()){
        bankStream = ifstream { TAAM_root + "/empty_TAAM_databank.txt" };
    }
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

    calculator.setAtomicFormfactorManager(formfactorCalculator);
}
