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
};
