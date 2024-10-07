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


// Forward declarations
void calculateSfTaamMinimal(
    const Crystal& crystal,
    const vector<Vector3i>& hkl,
    vector< complex<double> > &structureFactors);

void calculateSfIamMinimal(
    const Crystal& crystal,
    const vector<Vector3i>& hkl,
    vector< complex<double> >& structureFactors);
    

void DiscambWrapper::f_calc_hkl(const vector<Vector3i> &hkl, FCalcMethod method, vector<complex<double>> &sf){
    sf.clear();
    sf.resize(hkl.size());
    Crystal crystal;
    get_crystal(crystal);
    string table = get_discamb_table_string();
    switch (method)
    {
        case FCalcMethod::IAM: {
            auto calculator = get_IAM_calculator(crystal, table);
            calculator.setAnoumalous(mAnomalous);
            calculator.calculateStructureFactors(hkl, sf);
            break;
        }
        case FCalcMethod::TAAM: {
            auto calculator = get_TAAM_calculator(crystal);
            calculator.setAnoumalous(mAnomalous);
            calculator.calculateStructureFactors(hkl, sf);
            break;
        }
        default: {
            on_error::throwException("Invalid method requested");
            break;
        }
    }
}

vector<complex<double>> DiscambWrapper::f_calc(const double d_min, FCalcMethod method){
    vector<Vector3i> hkl;
    get_hkl(d_min, hkl);
    vector<complex<double>> sf;
    f_calc_hkl(hkl, method, sf);
    return sf;
}


void DiscambWrapper::get_crystal(Crystal &crystal){
    py::tuple cell_params = mStructure.attr("unit_cell")().attr("parameters")();
    crystal.unitCell.set(
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
    crystal.spaceGroup.set(symops);

    crystal.atoms.clear();
    mAnomalous.clear();
    for (auto scatterer_py : mStructure.attr("scatterers")()){
        crystal.atoms.push_back(AtomInCrystal());

        py::tuple xyz_py = scatterer_py.attr("site");
        crystal.atoms.back().coordinates[0] = xyz_py[0].cast<float>();
        crystal.atoms.back().coordinates[1] = xyz_py[1].cast<float>();
        crystal.atoms.back().coordinates[2] = xyz_py[2].cast<float>();

        crystal.atoms.back().coordinates_sigma[0] = 0.0;
        crystal.atoms.back().coordinates_sigma[1] = 0.0;
        crystal.atoms.back().coordinates_sigma[2] = 0.0;

        crystal.atoms.back().coordinates_precision[0] = 0.0;
        crystal.atoms.back().coordinates_precision[1] = 0.0;
        crystal.atoms.back().coordinates_precision[2] = 0.0;

        crystal.atoms.back().adp.clear();
        if (scatterer_py.attr("flags").attr("use_u_iso")().cast<bool>()){
            crystal.atoms.back().adp.push_back(scatterer_py.attr("u_iso").cast<float>()); // U11
        }
        else {
            py::tuple u_star = scatterer_py.attr("u_star");
            crystal.atoms.back().adp.push_back(u_star[0].cast<float>()); // U11
            crystal.atoms.back().adp.push_back(u_star[1].cast<float>()); // U22
            crystal.atoms.back().adp.push_back(u_star[2].cast<float>()); // U33
            crystal.atoms.back().adp.push_back(u_star[3].cast<float>()); // U12
            crystal.atoms.back().adp.push_back(u_star[4].cast<float>()); // U13
            crystal.atoms.back().adp.push_back(u_star[5].cast<float>()); // U23
        }
        crystal.atoms.back().adp_sigma.clear();
        crystal.atoms.back().adp_precision.clear();
        for(int i = 0; i < crystal.atoms.back().adp.size(); i++){
            crystal.atoms.back().adp_sigma.push_back(0.0);
            crystal.atoms.back().adp_precision.push_back(0.0);
        }
        
        crystal.atoms.back().label = scatterer_py.attr("scattering_type").cast<string>();

        crystal.atoms.back().occupancy = scatterer_py.attr("occupancy").cast<float>();
        crystal.atoms.back().multiplicity = scatterer_py.attr("multiplicity")().cast<float>();
        crystal.atoms.back().occupancy_sigma = 0.0;

        // This seems to be unused anyway
        crystal.atoms.back().siteSymetry.clear();
        crystal.atoms.back().siteSymetry.push_back(SpaceGroupOperation());

        mAnomalous.push_back(complex<double> {
            scatterer_py.attr("fp").cast<double>(),
            scatterer_py.attr("fdp").cast<double>()
        });
    }

    crystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;
    // Assume all adps use the same convention
    if (mStructure.attr("use_u_iso")().attr("__getitem__")(0).cast<bool>()){
        crystal.adpConvention = structural_parameters_convention::AdpConvention::U_cif;
    }
    else{
        crystal.adpConvention = structural_parameters_convention::AdpConvention::U_star;
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

void DiscambWrapper::update_atoms(Crystal &crystal){
    int idx = 0;
    for (auto scatterer_py : mStructure.attr("scatterers")()){

        py::tuple xyz_py = scatterer_py.attr("site");
        crystal.atoms[idx].coordinates[0] = xyz_py[0].cast<float>();
        crystal.atoms[idx].coordinates[1] = xyz_py[1].cast<float>();
        crystal.atoms[idx].coordinates[2] = xyz_py[2].cast<float>();

        crystal.atoms[idx].coordinates_sigma[0] = 0.0;
        crystal.atoms[idx].coordinates_sigma[1] = 0.0;
        crystal.atoms[idx].coordinates_sigma[2] = 0.0;

        crystal.atoms[idx].coordinates_precision[0] = 0.0;
        crystal.atoms[idx].coordinates_precision[1] = 0.0;
        crystal.atoms[idx].coordinates_precision[2] = 0.0;

        crystal.atoms[idx].adp.clear();
        crystal.atoms[idx].adp.push_back(scatterer_py.attr("u_iso").cast<float>()); // U11
        // crystal.atoms[idx].adp.push_back(scatterer_py.attr("u_iso").cast<float>()); // U22
        // crystal.atoms[idx].adp.push_back(scatterer_py.attr("u_iso").cast<float>()); // U33
        // crystal.atoms[idx].adp.push_back(0.0); // U12
        // crystal.atoms[idx].adp.push_back(0.0); // U13
        // crystal.atoms[idx].adp.push_back(0.0); // U23

        crystal.atoms[idx].adp_sigma.clear();
        crystal.atoms[idx].adp_sigma.push_back(0.0);
        // for(int i = 0; i < 6; i++){crystal.atoms[idx].adp_sigma.push_back(0.0);}

        crystal.atoms[idx].adp_precision.clear();
        crystal.atoms[idx].adp_precision.push_back(0.0);
        // for(int i = 0; i < 6; i++){crystal.atoms[idx].adp_precision.push_back(0.0);}
        
        crystal.atoms[idx].label = scatterer_py.attr("scattering_type").cast<string>();

        crystal.atoms[idx].occupancy = scatterer_py.attr("occupancy").cast<float>();
        crystal.atoms[idx].multiplicity = scatterer_py.attr("multiplicity")().cast<float>();
        crystal.atoms[idx].occupancy_sigma = 0.0;

        // This seems to be unused anyway
        crystal.atoms[idx].siteSymetry.clear();
        crystal.atoms[idx].siteSymetry.push_back(SpaceGroupOperation());

        mAnomalous[idx] = complex<double> {
            scatterer_py.attr("fp").cast<double>(),
            scatterer_py.attr("fdp").cast<double>()
        };
        idx++;
    }
}

string DiscambWrapper::get_discamb_table_string(){
    string cctbx_table_string = mStructure.attr("scattering_type_registry_params").attr("table").cast<string>();
    string discamb_table_string;
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
    else if (cctbx_table_string == "xray"){
        throw std::invalid_argument("Scattering table \"" + cctbx_table_string + "\" is not recognized.");
    }
    else if (cctbx_table_string == "neutron"){
        throw std::invalid_argument("Scattering table \"" + cctbx_table_string + "\" is not recognized.");
    }
    else if (cctbx_table_string == "n_gaussian"){
        throw std::invalid_argument("Scattering table \"" + cctbx_table_string + "\" is not recognized.");
    }
    return discamb_table_string;
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

    ifstream bankStream {"data/empty_TAAM_databank.txt"};
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
