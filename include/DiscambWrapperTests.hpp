#pragma once

#include "DiscambWrapper.hpp"

class DiscambWrapperTests : public DiscambWrapper{

    using DiscambWrapper::DiscambWrapper;

    public:
        void test_get_crystal(int n_iter);
        void test_update_atoms(int n_iter);
};
