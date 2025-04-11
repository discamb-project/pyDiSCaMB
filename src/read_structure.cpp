#include "read_structure.hpp"

#include "assert.hpp"
#include "scattering_table.hpp"

using namespace std;
using namespace discamb;

namespace py = pybind11;

Crystal crystal_from_xray_structure(const py::object &structure) {
    Crystal crystal;
    py::tuple cell_params = structure.attr("unit_cell")().attr("parameters")();
    crystal.unitCell.set(cell_params[0].cast<double>(),
                         cell_params[1].cast<double>(),
                         cell_params[2].cast<double>(),
                         cell_params[3].cast<double>(),
                         cell_params[4].cast<double>(),
                         cell_params[5].cast<double>());

    vector<SpaceGroupOperation> symops;
    for (auto symop_py :
         structure.attr("crystal_symmetry")().attr("space_group")()) {
        string symop_str = symop_py.attr("as_xyz")().cast<string>();
        SpaceGroupOperation symop{symop_str};
        symops.push_back(symop);
    }
    crystal.spaceGroup.set(symops);

    int num_atoms = structure.attr("scatterers")().attr("size")().cast<int>();
    crystal.atoms.clear();
    crystal.atoms.assign(num_atoms, AtomInCrystal());

    crystal.xyzCoordinateSystem =
        structural_parameters_convention::XyzCoordinateSystem::fractional;
    crystal.adpConvention =
        structural_parameters_convention::AdpConvention::U_star;

    // Set atoms
    update_crystal_from_xray_structure(crystal, structure);
    return crystal;
}

vector<complex<double>> anomalous_from_xray_structure(
    const py::object &structure) {
    vector<complex<double>> out;
    out.resize(structure.attr("scatterers")().attr("size")().cast<int>());
    update_anomalous_from_xray_structure(out, structure);
    return out;
}

string table_from_xray_structure(const py::object &structure) {
    auto py_table_string = structure.attr("get_scattering_table")();
    if (py::isinstance<py::none>(py_table_string)) {
        return table_alias("xray");
    }
    return table_alias(py_table_string.cast<string>());
}

void update_crystal_from_xray_structure(discamb::Crystal &crystal,
                                        const py::object &structure) {
    assert(crystal.atoms.size() ==
           structure.attr("scatterers")().attr("size")().cast<int>());

    int idx = 0;
    for (auto scatterer_py : structure.attr("scatterers")()) {
        // Type/label
        if (crystal.atoms[idx].label.size()) {
            assert(crystal.atoms[idx].label ==
                   scatterer_py.attr("label").cast<string>());
            assert(crystal.atoms[idx].type ==
                   scatterer_py.attr("scattering_type").cast<string>());
        } else {
            crystal.atoms[idx].type =
                scatterer_py.attr("scattering_type").cast<string>();
            crystal.atoms[idx].label =
                scatterer_py.attr("label").cast<string>();
        }

        // Site
        py::tuple xyz_py = scatterer_py.attr("site");
        crystal.atoms[idx].coordinates[0] = xyz_py[0].cast<float>();
        crystal.atoms[idx].coordinates[1] = xyz_py[1].cast<float>();
        crystal.atoms[idx].coordinates[2] = xyz_py[2].cast<float>();

        crystal.atoms[idx].coordinates_sigma[0] = 0.0;
        crystal.atoms[idx].coordinates_sigma[1] = 0.0;
        crystal.atoms[idx].coordinates_sigma[2] = 0.0;

        crystal.atoms[idx].coordinates_precision[0] = 0.0;
        crystal.atoms[idx].coordinates_precision[1] = 0.0;
        crystal.atoms[idx].coordinates_precision[2] = 0.0;

        // ADP
        crystal.atoms[idx].adp.clear();
        if (scatterer_py.attr("flags").attr("use_u_iso")().cast<bool>()) {
            crystal.atoms[idx].adp.push_back(
                scatterer_py.attr("u_iso").cast<float>());  // U11
        } else {
            py::tuple u_star = scatterer_py.attr("u_star");
            // TODO check if cctbx has the same order. Comments below are
            // discamb order
            crystal.atoms[idx].adp.push_back(u_star[0].cast<float>());  // U11
            crystal.atoms[idx].adp.push_back(u_star[1].cast<float>());  // U22
            crystal.atoms[idx].adp.push_back(u_star[2].cast<float>());  // U33
            crystal.atoms[idx].adp.push_back(u_star[3].cast<float>());  // U12
            crystal.atoms[idx].adp.push_back(u_star[4].cast<float>());  // U13
            crystal.atoms[idx].adp.push_back(u_star[5].cast<float>());  // U23
        }
        crystal.atoms[idx].adp_sigma.clear();
        crystal.atoms[idx].adp_precision.clear();
        for (int i = 0; i < crystal.atoms[idx].adp.size(); i++) {
            crystal.atoms[idx].adp_sigma.push_back(0.0);
            crystal.atoms[idx].adp_precision.push_back(0.0);
        }

        // occupancy
        crystal.atoms[idx].occupancy =
            scatterer_py.attr("occupancy").cast<float>();
        crystal.atoms[idx].multiplicity =
            scatterer_py.attr("multiplicity")().cast<float>();
        crystal.atoms[idx].occupancy_sigma = 0.0;

        // Site symmetry
        // This seems to be unused anyway
        crystal.atoms[idx].siteSymetry.clear();
        crystal.atoms[idx].siteSymetry.push_back(SpaceGroupOperation());

        idx++;
    }
}

void update_anomalous_from_xray_structure(vector<complex<double>> &anomalous,
                                          const py::object &structure) {
    assert(anomalous.size() ==
           structure.attr("scatterers")().attr("size")().cast<int>());
    int idx = 0;
    for (auto scatterer_py : structure.attr("scatterers")()) {
        anomalous[idx++] =
            complex<double>{scatterer_py.attr("fp").cast<double>(),
                            scatterer_py.attr("fdp").cast<double>()};
    }
}
