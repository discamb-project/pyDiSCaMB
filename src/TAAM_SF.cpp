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

        // TODO Basically just copy the above block for all the remaining attributes

    }
}

string test_func(const py::object structure){
    Crystal crystal;
    cctbx_model_to_discamb_crystal(structure, crystal);
    return "a: " + to_string(crystal.spaceGroup.nSymmetryOperations());
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
