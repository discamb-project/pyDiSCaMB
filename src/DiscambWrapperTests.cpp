#include "DiscambWrapperTests.hpp"

#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"

#include <map>

using namespace std;
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

vector<complex<double>> DiscambWrapperTests::test_f_calc(
    vector<string> atom_labels,
    vector<vector<double>> a,
    vector<vector<double>> b,
    double d_min
){
    vector<Vector3i> hkl;
    get_hkl(d_min, hkl);
    Crystal crystal;
    get_crystal(crystal);
    vector<complex<double>> sf;

    map<string, NGaussianFormFactor> scatterers;
    for (int i = 0; i < atom_labels.size(); i++){
        scatterers[atom_labels[i]] = NGaussianFormFactor(
            atom_labels[i],
            a[i],
            b[i],
            0.0
        );
    }

    shared_ptr<AtomicFormFactorCalculationsManager> formfactorCalculator;
    formfactorCalculator = shared_ptr<AtomicFormFactorCalculationsManager>(
        new IamFormFactorCalculationsManager(crystal, scatterers));
    AnyScattererStructureFactorCalculator calculator(crystal);
    calculator.setAtomicFormfactorManager(formfactorCalculator);

    calculator.calculateStructureFactors(hkl, sf);
    return sf;
}
