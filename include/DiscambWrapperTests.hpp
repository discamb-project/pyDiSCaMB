#pragma once

#include "DiscambWrapper.hpp"


class DiscambWrapperTests : public DiscambWrapper{

    using DiscambWrapper::DiscambWrapper;

    public:
        void test_get_crystal(int n_iter);
        void test_update_atoms(int n_iter);
        std::vector<std::complex<double>> test_f_calc(
            std::vector<std::string> atom_labels,
            std::vector<std::vector<double>> a,
            std::vector<std::vector<double>> b,
            double d_min
        );
        double get_f_calc_runtime(int n_iter, double d_min = 2.0);
        double get_f_calc_runtime_with_atom_updates(int n_iter, double d_min = 2.0);
};
