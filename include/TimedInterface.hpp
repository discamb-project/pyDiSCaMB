#pragma once

#include <string>
#include <vector>

#include "PythonInterface.hpp"

using Runtime = std::pair<std::string, double>;

class TimedInterface : public PythonInterface {
   public:
    TimedInterface(py::object& structure, py::dict params)
        : PythonInterface(structure, params) {
        init();
        update_structure(structure);
    }
    std::vector<Runtime> get_runtimes() const;

    // clang-format off
    void init() override;
    discamb::SfCalculator* init_calc(discamb::Crystal crystal, nlohmann::json calculator_parameters) override;
    void set_indices(py::object& indices) override;
    void set_d_min(const double d_min) override;
    void update_structure(py::object& structure) override;
    std::vector<std::complex<double>> f_calc() override;
    std::vector<FCalcDerivatives> d_f_calc_d_params() override;
    std::vector<discamb::TargetFunctionAtomicParamDerivatives> d_target_d_params(std::vector<std::complex<double>> d_target_d_f_calc) override;
    std::vector<discamb::TargetFunctionAtomicParamDerivatives> selected_d_target_d_params(std::vector<std::complex<double>> d_target_d_f_calc, bool site, bool adp, bool occupancy, bool fp) override;
    void update_calculator() override;
    // clang-format on

   private:
    std::vector<Runtime> mRuntimes;
};
