#include "DiscambWrapperTests.hpp"

#include "discamb/CrystalStructure/Crystal.h"

using namespace discamb;

void DiscambWrapperTests::test_get_crystal(int n_iter){
    for (int i = 0; i < n_iter; i++){
        Crystal crystal;
        get_crystal(crystal);
    }
}

void DiscambWrapperTests::test_update_atoms(int n_iter){
    Crystal crystal;
    get_crystal(crystal);
    for (int i = 0; i < n_iter; i++){
        update_atoms(crystal);
    }
}
