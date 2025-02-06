#include "tests.hpp"
#include "read_structure.hpp"

#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"

#include <map>

#include "assert.hpp"

using namespace std;
using namespace discamb;

namespace py = pybind11;

vector<complex<double>> f_calc_custom_gaussian_parameters(
    py::object structure,
    vector<vector<int>> hkl,
    vector<string> atom_labels,
    vector<vector<double>> a,
    vector<vector<double>> b,
    vector<double> c
){
    int i;
    assert(atom_labels.size() == a.size());
    assert(atom_labels.size() == b.size());
    assert(atom_labels.size() == c.size());
    for (i = 0; i < a.size(); i++){
        assert(a[i].size() == b[i].size());
    }
    for (i = 0; i < hkl.size(); i++){
        assert(hkl[i].size() == 3);
    }

    Crystal crystal = crystal_from_xray_structure(structure);

    int j;
    bool found;
    for (i = 0; i < crystal.atoms.size(); i++){
        found = false;
        for (j = 0; j < atom_labels.size(); j++){
            if (crystal.atoms[i].type == atom_labels[j]){
                found = true;
                break;
            }
        }
        if (!found) throw AssertionError("crystal.atoms[i].type in atom_labels", __FILE__, __LINE__);
    }

    vector<Vector3i> hkl_vector3i;
    for (i = 0; i < hkl.size(); i++){
        hkl_vector3i.push_back(Vector3i {hkl[i][0], hkl[i][1], hkl[i][2]});
    }

    vector<complex<double>> sf;
    sf.resize(hkl.size());

    map<string, NGaussianFormFactor> scatterers;
    for (i = 0; i < atom_labels.size(); i++){
        scatterers[atom_labels[i]] = NGaussianFormFactor(
            atom_labels[i],
            a[i],
            b[i],
            c[i]
        );
    }

    shared_ptr<AtomicFormFactorCalculationsManager> formfactorCalculator;
    formfactorCalculator = shared_ptr<AtomicFormFactorCalculationsManager>(
        new IamFormFactorCalculationsManager(crystal, scatterers));
    AnyScattererStructureFactorCalculator calculator(crystal);
    calculator.setAtomicFormfactorManager(formfactorCalculator);

    calculator.calculateStructureFactors(hkl_vector3i, sf);
    return sf;
}
