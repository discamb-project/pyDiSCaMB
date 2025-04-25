import pytest
from cctbx.eltbx.chemical_elements import proper_caps_list as elements
from pydiscamb import DiscambWrapper, FCalcMethod
from pydiscamb.taam_parameters import is_MATTS_installed


def _get_single_element_structure(element: str):
    from cctbx.development import random_structure as cctbx_random_structure
    from cctbx.sgtbx import space_group_info

    group = space_group_info(19)
    xrs = cctbx_random_structure.xray_structure(
        space_group_info=group,
        elements=[element],
        general_positions_only=False,
        use_u_iso=True,
        random_u_iso=False,
        random_occupancy=False,
    )
    xrs.scattering_type_registry(table="wk1995")
    return xrs


@pytest.mark.parametrize(
    "el",
    elements()[:36],
)
def test_assignment_light_element(el: str):
    xrs = _get_single_element_structure(el)

    w1 = DiscambWrapper(xrs, FCalcMethod.IAM)
    fc1 = w1.f_calc(2)
    w2 = DiscambWrapper(xrs, FCalcMethod.TAAM)
    fc2 = w2.f_calc(2)
    assert pytest.approx(fc1, rel=0.01) == fc2


@pytest.mark.parametrize(
    "el",
    elements()[36:86],
)
def test_assignment_heavy_element(el: str):
    xrs = _get_single_element_structure(el)

    w1 = DiscambWrapper(xrs, FCalcMethod.IAM)
    fc1 = w1.f_calc(2)
    w2 = DiscambWrapper(xrs, FCalcMethod.TAAM)
    fc2 = w2.f_calc(2)
    assert pytest.approx(fc1) == fc2


class TestAssignmentLogUnrecognizedStructure:

    def get_structure(self, atom_range: slice):
        from cctbx.development import random_structure as cctbx_random_structure
        from cctbx.sgtbx import space_group_info

        group = space_group_info(19)
        xrs = cctbx_random_structure.xray_structure(
            space_group_info=group,
            elements=elements()[atom_range],
            general_positions_only=False,
            use_u_iso=True,
            random_u_iso=False,
            random_occupancy=False,
        )
        xrs.scattering_type_registry(table="wk1995")
        return xrs

    def test_full(self, tmp_path):
        logfile = tmp_path / "assignment.log"

        xrs = self.get_structure(slice(0, 86, None))
        w = DiscambWrapper(xrs, FCalcMethod.TAAM, assignment_info=str(logfile))

        with logfile.open("r") as f:
            assert f.readline() == "Atom type assigned to 0 of 86.\n"
            # fmt: off
            assert f.readline() == "Atoms with unassigned atom types, represented with standard IAM :\n"
            for i in range(36, 86):
                assert f.readline() == f"{elements()[i]}{i + 1}\n".rjust(8)
            assert f.readline() == "Atoms with unassigned atom types, represented with multipole model based IAM :\n"
            # fmt: on
            for i in range(36):
                assert f.readline() == f"{elements()[i]}{i + 1}\n".rjust(8)
            f.readline()  # Empty line at the bottom
            assert f.readline() == ""  # Check for EOF

    def test_multipole(self, tmp_path):
        logfile = tmp_path / "assignment.log"

        xrs = self.get_structure(slice(0, 36, None))
        w = DiscambWrapper(xrs, FCalcMethod.TAAM, assignment_info=str(logfile))

        with logfile.open("r") as f:
            assert f.readline() == "Atom type assigned to 0 of 36.\n"
            # fmt: off
            assert f.readline() == "Atoms with unassigned atom types, represented with multipole model based IAM :\n"
            # fmt: on
            for i in range(36):
                assert f.readline() == f"{elements()[i]}{i + 1}\n".rjust(8)
            f.readline()  # Empty line at the bottom
            assert f.readline() == ""  # Check for EOF

    def test_iam(self, tmp_path):
        logfile = tmp_path / "assignment.log"

        xrs = self.get_structure(slice(36, 86, None))
        w = DiscambWrapper(xrs, FCalcMethod.TAAM, assignment_info=str(logfile))

        with logfile.open("r") as f:
            assert f.readline() == "Atom type assigned to 0 of 50.\n"
            # fmt: off
            assert f.readline() == "Atoms with unassigned atom types, represented with standard IAM :\n"
            # fmt: on
            for i in range(50):
                assert f.readline() == f"{elements()[i + 36]}{i + 1}\n".rjust(8)
            f.readline()  # Empty line at the bottom
            assert f.readline() == ""  # Check for EOF


@pytest.mark.skipif(not is_MATTS_installed(), reason="Must have MATTS installed")
def test_assignment_log_tyrosine(tyrosine, tmp_path):
    logfile = tmp_path / "assignment.log"

    w = DiscambWrapper(tyrosine, FCalcMethod.TAAM, assignment_info=str(logfile))

    with logfile.open("r") as f:
        assert f.readline() == "Atom type assigned to 24 of 24.\n"
        # fmt: off
        assert f.readline() == "Atoms with assigned atom types and local coordinate systems:\n"
        assert f.readline() == '   pdb=" N   TYR A   4 "   N401a    Z pdb=" CA  TYR A   4 " X pdb=" H2  TYR A   4 "\n'
        assert f.readline() == '   pdb=" CA  TYR A   4 "    C414    X pdb=" N   TYR A   4 " Y pdb=" C   TYR A   4 "\n'
        assert f.readline() == '   pdb=" C   TYR A   4 "    C301    Z pdb=" CA  TYR A   4 " X pdb=" OXT TYR A   4 "\n'
        assert f.readline() == '   pdb=" O   TYR A   4 "    O101    X pdb=" C   TYR A   4 " Y pdb=" OXT TYR A   4 "\n'
        assert f.readline() == '   pdb=" CB  TYR A   4 "    C405    X pdb=" CG  TYR A   4 " Y pdb=" CA  TYR A   4 "\n'
        assert f.readline() == '   pdb=" CG  TYR A   4 "    C332    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CD2 TYR A   4 "\n'
        assert f.readline() == '   pdb=" CD1 TYR A   4 "    C330    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CE1 TYR A   4 "\n'
        assert f.readline() == '   pdb=" CD2 TYR A   4 "    C330    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CE2 TYR A   4 "\n'
        assert f.readline() == '   pdb=" CE1 TYR A   4 "    C330    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CZ  TYR A   4 "\n'
        assert f.readline() == '   pdb=" CE2 TYR A   4 "    C330    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CZ  TYR A   4 "\n'
        assert f.readline() == '   pdb=" CZ  TYR A   4 "   C341a    X average_position(pdb=" CG  TYR A   4 ",pdb=" CD1 TYR A   4 ",pdb=" CE1 TYR A   4 ",pdb=" CZ  TYR A   4 ",pdb=" CE2 TYR A   4 ",pdb=" CD2 TYR A   4 ") Y pdb=" CE2 TYR A   4 "\n'
        assert f.readline() == '   pdb=" OH  TYR A   4 "    O204    X pdb=" CZ  TYR A   4 " Y pdb=" HH  TYR A   4 "\n'
        assert f.readline() == '   pdb=" OXT TYR A   4 "    O101    X pdb=" C   TYR A   4 " Y pdb=" O   TYR A   4 "\n'
        assert f.readline() == '   pdb=" H1  TYR A   4 "    H105    Z pdb=" N   TYR A   4 " X pdb=" CA  TYR A   4 "\n'
        assert f.readline() == '   pdb=" H2  TYR A   4 "    H105    Z pdb=" N   TYR A   4 " X pdb=" CA  TYR A   4 "\n'
        assert f.readline() == '   pdb=" H3  TYR A   4 "    H105    Z pdb=" N   TYR A   4 " X pdb=" CA  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HA  TYR A   4 "    H103    Z pdb=" CA  TYR A   4 " X pdb=" N   TYR A   4 "\n'
        assert f.readline() == '   pdb=" HB2 TYR A   4 "    H102    Z pdb=" CB  TYR A   4 " X pdb=" CG  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HB3 TYR A   4 "    H102    Z pdb=" CB  TYR A   4 " X pdb=" CG  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HD1 TYR A   4 "    H104    Z pdb=" CD1 TYR A   4 " X pdb=" CE1 TYR A   4 "\n'
        assert f.readline() == '   pdb=" HD2 TYR A   4 "    H104    Z pdb=" CD2 TYR A   4 " X pdb=" CE2 TYR A   4 "\n'
        assert f.readline() == '   pdb=" HE1 TYR A   4 "    H104    Z pdb=" CE1 TYR A   4 " X pdb=" CZ  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HE2 TYR A   4 "    H104    Z pdb=" CE2 TYR A   4 " X pdb=" CZ  TYR A   4 "\n'
        assert f.readline() == '   pdb=" HH  TYR A   4 "    H114    Z pdb=" OH  TYR A   4 " X pdb=" CZ  TYR A   4 "\n'
        # fmt: on
        f.readline()  # Empty line at the bottom
        assert f.readline() == ""  # Check for EOF

@pytest.mark.skipif(not is_MATTS_installed(), reason="Must have MATTS installed")
def test_assignment_symmetry_wrapping(tmp_path):
    # Download a graphene cif
    import iotbx.cif
    import requests

    data = requests.get(
        "https://legacy.materialsproject.org/materials/mp-48/cif?type=symmetrized&download=true"
    )
    cif_str = data.content.decode("utf-8")

    cif = iotbx.cif.reader(input_string=cif_str)
    xrs = cif.build_crystal_structures()["C"]
    xrs.scattering_type_registry(table="it1992")

    logfile = tmp_path / "assignment.log"
    w = DiscambWrapper(xrs, FCalcMethod.TAAM, assignment_info=str(logfile))
    with logfile.open("r") as f:
        assert f.readline() == "Atom type assigned to 2 of 2.\n"


@pytest.mark.skipif(not is_MATTS_installed(), reason="Must have MATTS installed")
def test_assignment_break_when_changing_to_heavy_atom(tmp_path, ethane):
    # Ethene in a big P1 box.
    # All H should be H101, both C should be C401

    logfile = tmp_path / "assignment.log"
    w = DiscambWrapper(ethane, FCalcMethod.TAAM, assignment_info=str(logfile))
    with logfile.open("r") as f:
        assert f.readline() == "Atom type assigned to 8 of 8.\n"
        f.readline()
        assert "C1   C401" in f.readline()
        assert "C2   C401" in f.readline()
        assert "H1   H101" in f.readline()
        assert "H2   H101" in f.readline()
        assert "H3   H101" in f.readline()
        assert "H4   H101" in f.readline()
        assert "H5   H101" in f.readline()
        assert "H6   H101" in f.readline()

    # Now: change one of the carbons to gold.
    # The three hydrogens on the carbon should still be H101,
    # but the other three hydrogens should be spherical
    xrs = ethane.deep_copy_scatterers()
    xrs.scatterers()[1].label = "Au1"
    xrs.scatterers()[1].scattering_type = "Au"
    w = DiscambWrapper(xrs, FCalcMethod.TAAM, assignment_info=str(logfile))
    with logfile.open("r") as f:
        assert f.readline() == "Atom type assigned to 3 of 8.\n"
        f.readline()
        assert f.readline() == "   Au1\n"
        f.readline()
        assert f.readline() == "    C1\n"
        assert f.readline() == "    H4\n"
        assert f.readline() == "    H5\n"
        assert f.readline() == "    H6\n"
        f.readline()
        assert f.readline() == "    H1        H101    Z C1 X Au1\n"
        assert f.readline() == "    H2        H101    Z C1 X Au1\n"
        assert f.readline() == "    H3        H101    Z C1 X Au1\n"


class TestAssignmentDict:
    def test_spherical(self):
        xrs = _get_single_element_structure("Na")
        w = DiscambWrapper(xrs, FCalcMethod.TAAM)
        assert hasattr(w, "atom_type_assignment")
        assert isinstance(w.atom_type_assignment, dict)
        assert len(w.atom_type_assignment) == 1
        assert w.atom_type_assignment["Na1"] == ("Na000", "")

    def test_untyped(self):
        xrs = _get_single_element_structure("Au")
        w = DiscambWrapper(xrs, FCalcMethod.TAAM)
        assert hasattr(w, "atom_type_assignment")
        assert isinstance(w.atom_type_assignment, dict)
        assert len(w.atom_type_assignment) == 1
        assert w.atom_type_assignment["Au1"] == ("", "")

    def test_multipolar(self, ethane):
        w = DiscambWrapper(ethane, FCalcMethod.TAAM)
        assert hasattr(w, "atom_type_assignment")
        assert isinstance(w.atom_type_assignment, dict)
        assert len(w.atom_type_assignment) == ethane.scatterers().size()
        assert all(a for a, b in w.atom_type_assignment.values())
        assert all(b for a, b in w.atom_type_assignment.values())

    def test_IAM(self, ethane):
        w = DiscambWrapper(ethane, FCalcMethod.IAM)
        assert hasattr(w, "atom_type_assignment")
        assert isinstance(w.atom_type_assignment, dict)
        assert len(w.atom_type_assignment) == ethane.scatterers().size()
        assert not any(a for a, b in w.atom_type_assignment.values())
        assert not any(b for a, b in w.atom_type_assignment.values())

    def test_explicit_file_exists(self, random_structure, tmp_path):
        csv = tmp_path / "assignment_log.csv"
        w = DiscambWrapper(random_structure, FCalcMethod.TAAM, assignment_csv=str(csv))
        assert csv.exists()
        assert hasattr(w, "atom_type_assignment")
        assert isinstance(w.atom_type_assignment, dict)
        assert len(w.atom_type_assignment) == random_structure.scatterers().size()

    def test_csv_is_deleted(self, random_structure):
        import os, tempfile

        n_tmpfiles_before = len(os.listdir(tempfile.gettempdir()))
        w = DiscambWrapper(random_structure, FCalcMethod.TAAM)
        n_tmpfiles_after = len(os.listdir(tempfile.gettempdir()))
        assert n_tmpfiles_after == n_tmpfiles_before
        assert hasattr(w, "atom_type_assignment")
        assert isinstance(w.atom_type_assignment, dict)
        assert len(w.atom_type_assignment) == random_structure.scatterers().size()
