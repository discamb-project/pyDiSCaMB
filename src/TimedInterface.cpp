#include "TimedInterface.hpp"

#include <chrono>

using namespace std;
using namespace discamb;
using namespace chrono;

namespace py = pybind11;

std::vector<Runtime> TimedInterface::get_runtimes() const { return mRuntimes; }

// clang-format off
void TimedInterface::init() {
    auto runtime_start = high_resolution_clock::now();
    PythonInterface::init();
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("init"), delta.count()});
}

discamb::SfCalculator* TimedInterface::init_calc(discamb::Crystal crystal, nlohmann::json calculator_parameters) {
    auto runtime_start = high_resolution_clock::now();
    discamb::SfCalculator* res = PythonInterface::init_calc(crystal, calculator_parameters);
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("init_calc"), delta.count()});
    return res;
}

void TimedInterface::set_indices(py::object& indices) {
    auto runtime_start = high_resolution_clock::now();
    PythonInterface::set_indices(indices);
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("set_indices"), delta.count()});
}

void TimedInterface::set_d_min(const double d_min) {
    auto runtime_start = high_resolution_clock::now();
    PythonInterface::set_d_min(d_min);
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("set_d_min"), delta.count()});
}

void TimedInterface::update_structure(py::object& structure) {
    auto runtime_start = high_resolution_clock::now();
    PythonInterface::update_structure(structure);
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("update_structure"), delta.count()});
}

std::vector<std::complex<double>> TimedInterface::f_calc() {
    auto runtime_start = high_resolution_clock::now();
    std::vector<std::complex<double>> res = PythonInterface::f_calc();
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("f_calc"), delta.count()});
    return res;
}

std::vector<FCalcDerivatives> TimedInterface::d_f_calc_d_params() {
    auto runtime_start = high_resolution_clock::now();
    std::vector<FCalcDerivatives> res = PythonInterface::d_f_calc_d_params();
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("d_f_calc_d_params"), delta.count()});
    return res;
}

std::vector<discamb::TargetFunctionAtomicParamDerivatives> TimedInterface::d_target_d_params(std::vector<std::complex<double>> d_target_d_f_calc) {
    auto runtime_start = high_resolution_clock::now();
    std::vector<discamb::TargetFunctionAtomicParamDerivatives> res = PythonInterface::d_target_d_params(d_target_d_f_calc);
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("d_target_d_params"), delta.count()});
    return res;
}

std::vector<discamb::TargetFunctionAtomicParamDerivatives> TimedInterface::selected_d_target_d_params(std::vector<std::complex<double>> d_target_d_f_calc, bool site, bool adp, bool occupancy, bool fp) {
    auto runtime_start = high_resolution_clock::now();
    std::vector<discamb::TargetFunctionAtomicParamDerivatives> res = PythonInterface::selected_d_target_d_params(d_target_d_f_calc, site, adp, occupancy, fp);
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("selected_d_target_d_params"), delta.count()});
    return res;
}

void TimedInterface::update_calculator() {
    auto runtime_start = high_resolution_clock::now();
    PythonInterface::update_calculator();
    auto runtime_end = high_resolution_clock::now();
    duration<double> delta = runtime_end - runtime_start;
    mRuntimes.push_back({std::string("update_calculator"), delta.count()});
}
// clang-format on
