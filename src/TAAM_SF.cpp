#include "pybind11/pybind11.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/UnitCell.h"
#include "discamb/CrystalStructure/SpaceGroup.h"
#include "discamb/CrystalStructure/SpaceGroupOperation.h"

#include "discamb_wrapper.hpp"


namespace py = pybind11;

using namespace std;
using namespace discamb;


void cctbx_model_to_discamb_crystal(const py::object model, Crystal &crystal){
    py::object structure = model.attr("get_xray_structure")();
    py::tuple cell_params = structure.attr("unit_cell")().attr("parameters")();
    UnitCell cell {
        cell_params[0].cast<double>(), 
        cell_params[1].cast<double>(), 
        cell_params[2].cast<double>(), 
        cell_params[3].cast<double>(), 
        cell_params[4].cast<double>(), 
        cell_params[5].cast<double>()
        };
    crystal.unitCell = cell;

    vector<SpaceGroupOperation> symops;
    for (auto symop_py : structure.attr("crystal_symmetry")().attr("space_group")()){
        string symop_str = symop_py.attr("as_xyz")().cast<string>();
        SpaceGroupOperation symop {symop_str};
        symops.push_back(symop);
    }
    SpaceGroup space_group;
    space_group.set(symops);
    crystal.spaceGroup = space_group;

    vector<AtomInCrystal> atoms;
    for (auto atom_py : model.attr("get_atoms")()){
        AtomInCrystal atom;

        py::tuple xyz_py = atom_py.attr("xyz");
        atom.coordinates[0] = xyz_py[0].cast<double>();
        atom.coordinates[1] = xyz_py[1].cast<double>();
        atom.coordinates[2] = xyz_py[2].cast<double>();

        py::tuple sigxyz_py = atom_py.attr("sigxyz");
        atom.coordinates_sigma[0] = sigxyz_py[0].cast<double>();
        atom.coordinates_sigma[1] = sigxyz_py[1].cast<double>();
        atom.coordinates_sigma[2] = sigxyz_py[2].cast<double>();

        // assume single adp for now
        atom.adp.resize(1);
        atom.adp.push_back(atom_py.attr("b").cast<double>());
        atom.adp_sigma.resize(1);
        atom.adp_sigma.push_back(atom_py.attr("sigb").cast<double>());
        
        atom.label = atom_py.attr("name").cast<string>();

        atom.occupancy = atom_py.attr("occ").cast<double>();
        atom.occupancy_sigma = atom_py.attr("sigocc").cast<double>();

        // This seems to be unused anyway
        atom.siteSymetry.resize(1);
        SpaceGroupOperation identity = SpaceGroupOperation();
        atom.siteSymetry.push_back(identity);

        atoms.push_back(atom);
    }
    crystal.atoms = atoms;
    
    // Not sure if the coordinate format can change in the python object?
    crystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::cartesian;
}

void cctbx_model_to_hkl(py::object model, double d, vector<Vector3i> &hkl){
    // assume anomolous_flag is False
    py::object miller_py = model.attr("get_xray_structure")().attr("build_miller_set")(false, d);

    hkl.resize(0);
    for (auto hkl_py_auto : miller_py.attr("indices")().attr("as_vec3_double")()){
        py::tuple hkl_py = hkl_py_auto.cast<py::tuple>();
        Vector3i hkl_i {
            static_cast<int>(hkl_py[0].cast<double>()),
            static_cast<int>(hkl_py[1].cast<double>()),
            static_cast<int>(hkl_py[2].cast<double>())
        };
        hkl.push_back(hkl_i);
    }
}

string test_func(const py::object model){
    Crystal crystal;
    vector<Vector3i> hkl;
    cctbx_model_to_discamb_crystal(model, crystal);
    cctbx_model_to_hkl(model, 0.15, hkl);

    vector< complex<double> > structureFactors;
    structureFactors.reserve(hkl.size());
    calculateSfTaamMinimal(crystal, hkl, structureFactors);

    int ind = 10;
    complex<double> sf = structureFactors[ind];
    double sf_real = static_cast<double>(sf.real());
    double sf_imag = static_cast<double>(sf.imag());
    return "hkl: " + to_string(hkl[ind][0]) + to_string(hkl[ind][1]) + to_string(hkl[ind][2]) + "\nsf: " + to_string(sf_real) + " + " + to_string(sf_imag) + "j";
}

PYBIND11_MODULE(_taam_sf, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: taam_sf

        .. autosummary::
           :toctree: _generate

           get_discamb_version
    )pbdoc";

    m.def("get_discamb_version", &get_discamb_version, R"pbdoc(
        Get the version string for DiSCaMB
    )pbdoc");

    m.def("test", &test_func, "placeholder");
}
