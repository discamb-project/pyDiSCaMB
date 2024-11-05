#include "DiscambWrapperTests.hpp"

#include "discamb/BasicUtilities/Timer.h"

#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"

#include <map>

using namespace std;
using namespace discamb;

void DiscambWrapperTests::test_get_crystal(int n_iter){
    for (int i = 0; i < n_iter; i++){
        init_crystal();
    }
}

void DiscambWrapperTests::test_update_atoms(int n_iter){
    for (int i = 0; i < n_iter; i++){
        update_atoms();
    }
}

vector<complex<double>> DiscambWrapperTests::test_f_calc(
    vector<string> atom_labels,
    vector<vector<double>> a,
    vector<vector<double>> b,
    double d_min
){
    set_d_min(d_min);

    vector<complex<double>> sf;
    sf.resize(mHkl.size());

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
        new IamFormFactorCalculationsManager(mCrystal, scatterers));
    AnyScattererStructureFactorCalculator calculator(mCrystal);
    calculator.setAtomicFormfactorManager(formfactorCalculator);

    calculator.calculateStructureFactors(mHkl, sf);
    return sf;
}

double DiscambWrapperTests::get_f_calc_runtime(int n_iter, double d_min){
    CpuTimer timer;
    
    set_d_min(d_min);
    vector<complex<double>> sf;
    sf.resize(mHkl.size());
    update();
    timer.start();
    for (int i = 0; i < n_iter; i++){
        mCalculator.calculateStructureFactors(mHkl, sf);
    }
    return timer.stop();
}
double DiscambWrapperTests::get_f_calc_runtime_with_atom_updates(int n_iter, double d_min){
    CpuTimer timer;

    set_d_min(d_min);
    vector<complex<double>> sf;
    sf.resize(mHkl.size());
    timer.start();
    for (int i = 0; i < n_iter; i++){
        update();
        mCalculator.calculateStructureFactors(mHkl, sf);
    }
    return timer.stop();
}
