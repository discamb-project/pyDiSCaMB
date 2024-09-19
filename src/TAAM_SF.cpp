#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/UnitCell.h"
#include "discamb/CrystalStructure/SpaceGroup.h"
#include "discamb/CrystalStructure/SpaceGroupOperation.h"

#include "discamb/BasicUtilities/discamb_version.h"

#include "discamb_wrapper.hpp"


namespace py = pybind11;

using namespace std;
using namespace discamb;


void cctbx_model_to_discamb_crystal(const py::object structure, Crystal &crystal){
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
    for (auto scatterer_py : structure.attr("scatterers")()){
        AtomInCrystal atom;

        py::tuple xyz_py = scatterer_py.attr("site");
        atom.coordinates[0] = xyz_py[0].cast<double>();
        atom.coordinates[1] = xyz_py[1].cast<double>();
        atom.coordinates[2] = xyz_py[2].cast<double>();

        atom.coordinates_sigma[0] = 0.0;
        atom.coordinates_sigma[1] = 0.0;
        atom.coordinates_sigma[2] = 0.0;

        // assume single adp for now
        atom.adp.resize(1);
        atom.adp.push_back(scatterer_py.attr("u_iso").cast<double>());
        atom.adp_sigma.resize(1);
        atom.adp_sigma.push_back(0.0);
        
        atom.label = scatterer_py.attr("scattering_type").cast<string>();

        atom.occupancy = scatterer_py.attr("occupancy").cast<double>();
        atom.multiplicity = scatterer_py.attr("multiplicity")().cast<double>();
        atom.occupancy_sigma = 0.0;

        // This seems to be unused anyway
        atom.siteSymetry.resize(1);
        SpaceGroupOperation identity = SpaceGroupOperation();
        atom.siteSymetry.push_back(identity);

        atoms.push_back(atom);
    }
    crystal.atoms = atoms;

    crystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;
    crystal.adpConvention = structural_parameters_convention::AdpConvention::U_cif;
}

void cctbx_model_to_hkl(py::object structure, double d, vector<Vector3i> &hkl){
    // assume anomolous_flag is False
    py::object miller_py = structure.attr("build_miller_set")(false, d);

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

vector<complex<double>> test_TAAM(const py::object structure, double d){
    Crystal crystal;
    vector<Vector3i> hkl;
    cctbx_model_to_discamb_crystal(structure, crystal);
    cctbx_model_to_hkl(structure, d, hkl);

    vector< complex<double> > structureFactors;
    structureFactors.reserve(hkl.size());
    calculateSfTaamMinimal(crystal, hkl, structureFactors);

    return structureFactors;
}

vector<complex<double>> test_IAM(const py::object structure, double d){
    Crystal crystal;
    vector<Vector3i> hkl;
    cctbx_model_to_discamb_crystal(structure, crystal);
    cctbx_model_to_hkl(structure, d, hkl);

    vector< complex<double> > structureFactors;
    structureFactors.reserve(hkl.size());
    calculateSfIamMinimal(crystal, hkl, structureFactors);

    return structureFactors;
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

    m.def("get_discamb_version", &discamb_version::version, R"pbdoc(
        Get the version string for DiSCaMB
    )pbdoc");

    m.def("test_TAAM", &test_TAAM, "placeholder docstring");
    m.def("test_IAM", &test_IAM, "placeholder docstring");
}
