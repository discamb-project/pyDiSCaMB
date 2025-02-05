#include "python_wrapper.hpp"
#include "scattering_table.hpp"

#include "discamb/MathUtilities/Vector3.h"

#include <utility>

#include "read_structure.hpp"

using namespace std;
using namespace discamb;

namespace py = pybind11;

SfCalculator *get_calculator(py::object structure, FCalcMethod method){
    Crystal crystal = crystal_from_xray_structure(structure);
    nlohmann::json calculator_params;
    switch (method)
    {
    case FCalcMethod::IAM: {
        calculator_params = {
            {"model", "iam"},
            {"electron scattering", false},
            {"table", table_alias(structure.attr("get_scattering_table")().cast<string>())},
        };
        break;
    }
    case FCalcMethod::TAAM: {
        calculator_params = {
            {"model", "matts"},
            {"electron scattering", table_alias(structure.attr("get_scattering_table")().cast<string>()).find("electron") != string::npos},
            {"bank path", py::module::import("pydiscamb.taam_parameters").attr("get_default_databank")().cast<string>()},
        };
        break;
    }
    default:
        break;
    }
    return SfCalculator::create(crystal, calculator_params);
}

DiscambWrapper::DiscambWrapper(py::object structure, FCalcMethod method) :
    mStructure(std::move(structure)),
    mDiscambCalculator(
        get_calculator(mStructure, method),
        crystal_from_xray_structure(mStructure),
        anomalous_from_xray_structure(mStructure)
    ) 
    {};

DiscambWrapper DiscambWrapper::from_TAAM_parameters(
    py::object structure,
    bool convert_to_electron_scattering,
    std::string bank_filepath,
    std::string assignment_log_filepath,
    std::string parameter_log_filepath,
    std::string multipolar_cif_output_filepath,
    double unit_cell_charge,
    bool perform_parameter_scaling_from_unit_cell_charge
){
    DiscambWrapper out = DiscambWrapper(structure, FCalcMethod::IAM);
    Crystal crystal = crystal_from_xray_structure(mStructure);
    vector<complex<double>> anomalous = anomalous_from_xray_structure(mStructure);
    nlohmann::json params {
        {"model", "matts"},
        {"electron scattering", convert_to_electron_scattering},
        {"bank path", bank_filepath},
        {"assignment info", assignment_log_filepath},
        {"parameters info", parameter_log_filepath},
        {"multipole cif", multipolar_cif_output_filepath},
        {"unit cell charge", unit_cell_charge},
        {"scale", perform_parameter_scaling_from_unit_cell_charge}
    };
    out.mDiscambCalculator = DiscambStructureFactorCalculator(SfCalculator::create(crystal, params), crystal, anomalous);
    return out;
}

void DiscambWrapper::set_indices(py::object indices){
    mDiscambCalculator.hkl.clear();
    for (auto hkl_py_auto : indices){
        py::tuple hkl_py = hkl_py_auto.cast<py::tuple>();
        mDiscambCalculator.hkl.push_back(Vector3i {
            hkl_py[0].cast<int>(),
            hkl_py[1].cast<int>(),
            hkl_py[2].cast<int>()
        });
    }
}

void DiscambWrapper::set_d_min(const double d_min){
    bool anomalous_flag = mStructure.attr("scatterers")().attr("count_anomalous")().cast<int>() != 0;
    py::object miller_py = mStructure.attr("build_miller_set")(anomalous_flag, d_min);

    set_indices(miller_py.attr("indices")());
}

vector<complex<double>> DiscambWrapper::f_calc(){
    return mDiscambCalculator.f_calc();
}
vector<complex<double>> DiscambWrapper::f_calc(const double d_min){
    set_d_min(d_min);
    return f_calc();
}

vector<FCalcDerivatives> DiscambWrapper::d_f_calc_d_params(){
    return mDiscambCalculator.d_f_calc_d_params();
}

FCalcDerivatives DiscambWrapper::d_f_calc_hkl_d_params(py::tuple hkl){
    return d_f_calc_hkl_d_params(hkl[0].cast<int>(), hkl[1].cast<int>(), hkl[2].cast<int>());
}

FCalcDerivatives DiscambWrapper::d_f_calc_hkl_d_params(int h, int k, int l){
    return mDiscambCalculator.d_f_calc_hkl_d_params(h, k, l);
}

vector<TargetFunctionAtomicParamDerivatives> DiscambWrapper::d_target_d_params(vector<complex<double>> d_target_d_f_calc){
    return mDiscambCalculator.d_target_d_params(d_target_d_f_calc);
}


vector<complex<double>> calculate_structure_factors(py::object structure, double d, FCalcMethod method){
    DiscambWrapper w {structure, method};
    vector< complex<double> > structureFactors = w.f_calc(d);
    return structureFactors;
}

vector<complex<double>> calculate_structure_factors_TAAM(py::object structure, double d){
    return calculate_structure_factors(structure, d, FCalcMethod::TAAM);
}

vector<complex<double>> calculate_structure_factors_IAM(py::object structure, double d){
    return calculate_structure_factors(structure, d, FCalcMethod::IAM);
}
