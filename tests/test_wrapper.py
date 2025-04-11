import pytest
from cctbx import miller
from cctbx.array_family import flex
from pydiscamb import DiscambWrapper, FCalcMethod
from pydiscamb.discamb_wrapper import DiscambWrapperCached



class TestInit:
    def test_simple(self, random_structure):
        w = DiscambWrapper(random_structure)

    @pytest.mark.parametrize("taam", [True, False])
    def test_default_init_with_kwargs(self, tyrosine, taam):
        m = FCalcMethod.TAAM if taam else FCalcMethod.IAM
        w1 = DiscambWrapper(tyrosine, m)
        w2 = DiscambWrapper(tyrosine, m, random_unused_kwarg=42)

        fc1 = w1.f_calc(2.0)
        fc2 = w2.f_calc(2.0)
        assert pytest.approx(list(fc1)) == list(fc2)

    @pytest.mark.slow
    @pytest.mark.parametrize(
            "method", 
            [
                # Multiple runs to ensure cached code does not interfere either
                FCalcMethod.TAAM,
                FCalcMethod.IAM, 
                FCalcMethod.IAM, 
                FCalcMethod.TAAM,
                FCalcMethod.IAM, 
                FCalcMethod.TAAM,
                ]
                )
    def test_no_memory_leak(self, tyrosine, method):
        import psutil
        
        def instatiate_wrapper(n):
            hkl = flex.miller_index(1000, (1, 2, 3))
            for i in range(n):
                w = DiscambWrapper(tyrosine, method)
                w.set_indices(hkl)
                assert len(w.hkl) == 1000
        instatiate_wrapper(1)
        
        mem_before = psutil.Process().memory_info()
        instatiate_wrapper(200)
        mem_after = psutil.Process().memory_info()
        assert mem_before.vms == mem_after.vms
        assert mem_before.rss == mem_after.rss

class TestFCalc:
    def test_simple(self, random_structure):
        w = DiscambWrapper(random_structure)
        sf = w.f_calc(4.0)
        assert isinstance(sf[0], complex)

    def test_d_min(self, random_structure):
        w = DiscambWrapper(random_structure)
        w.set_d_min(4.0)
        sf = w.f_calc()
        assert isinstance(sf[0], complex)

    def test_indices(self, random_structure):
        inds = [(0, 1, 2), (0, 1, 2), (10, 20, 30)]

        w = DiscambWrapper(random_structure)
        w.set_indices(inds)
        sf = w.f_calc()
        assert len(sf) == len(inds)
        assert isinstance(sf[0], complex)
        assert pytest.approx(sf[0]) == sf[1]
        assert pytest.approx(sf[0]) != sf[2]

    def test_no_indices(self, random_structure):
        w = DiscambWrapper(random_structure)
        assert w.hkl == []
        sf = w.f_calc()
        assert len(sf) == 0

    @pytest.mark.parametrize(["p", "dp"], [(False, True), (True, False), (True, True)])
    def test_anomalous_scattering(self, p: bool, dp: bool, random_structure):
        wrapper = DiscambWrapper(random_structure)
        sf_before = wrapper.f_calc(5)

        sf_before = random_structure.structure_factors(d_min=5).f_calc().data()
        rb = [sf.real for sf in sf_before]
        ib = [sf.imag for sf in sf_before]
        b = [abs(sf) for sf in sf_before]

        if p:
            random_structure.shake_fps()
        if dp:
            random_structure.shake_fdps()

        wrapper = DiscambWrapper(random_structure)
        sf_after = wrapper.f_calc(5)
        sf_after = random_structure.structure_factors(d_min=5).f_calc().data()
        ra = [sf.real for sf in sf_after]
        ia = [sf.imag for sf in sf_after]
        a = [abs(sf) for sf in sf_after]

        assert not pytest.approx(rb) == ra
        assert not pytest.approx(ib) == ia
        assert not pytest.approx(b) == a

    def test_types(self, tyrosine):
        fc_c = tyrosine.structure_factors(d_min=2).f_calc()

        w = DiscambWrapper(tyrosine)
        w.set_indices(fc_c.indices())
        assert isinstance(w.f_calc(), flex.complex_double)

        w = DiscambWrapper(tyrosine)
        assert isinstance(w.f_calc(fc_c), miller.array)

        with pytest.raises(
            ValueError, match="`miller_set` must be of type `cctbx.miller.set"
        ):
            w.f_calc("incorrect input")


class TestUpdateStructure:
    def test_simple(self, large_random_structure):
        wrapper = DiscambWrapper(large_random_structure)
        sf_before = wrapper.f_calc(5)
        large_random_structure.shake_sites_in_place(0.1)
        wrapper.update_structure(large_random_structure)
        sf_after = wrapper.f_calc(5)
        assert pytest.approx(list(sf_before)) != list(sf_after)

    def test_raise_on_different_atoms(self, random_structure, tyrosine):
        w = DiscambWrapper(random_structure)
        with pytest.raises(ValueError, match="Incompatible structures"):
            w.update_structure(tyrosine)

    def test_raise_on_different_unit_cell(
        self, random_structure, random_structure_u_iso
    ):
        assert "".join(s.label for s in random_structure.scatterers()) == "".join(
            s.label for s in random_structure_u_iso.scatterers()
        )
        assert random_structure_u_iso.unit_cell() != random_structure.unit_cell()

        w = DiscambWrapper(random_structure)
        with pytest.raises(ValueError, match="Incompatible structures"):
            w.update_structure(random_structure_u_iso)


class TestCache:
    def test_simple(self, random_structure):
        DiscambWrapperCached.__cache = {}
        # Show that, even when objects are created in an unavailable scope,
        # the object is still cached

        def init_wrapper():
            w = DiscambWrapperCached(random_structure)
            res = id(w)
            del w
            return res

        with pytest.raises(NameError):
            # Wrapper is not in this scope
            print(w)

        wrapper_ids = [
            init_wrapper(),
            init_wrapper(),
            init_wrapper(),
            init_wrapper(),
        ]
        assert all(wi == init_wrapper() for wi in wrapper_ids)

    @pytest.mark.xfail(
        reason="Assignment is very fast now. Update test to larger structure"
    )
    def test_timesave(self, lysozyme):
        DiscambWrapperCached.__cache = {}
        from time import perf_counter

        start = perf_counter()
        w1 = DiscambWrapperCached(lysozyme, FCalcMethod.TAAM)
        w1_time = perf_counter() - start

        start = perf_counter()
        w2 = DiscambWrapperCached(lysozyme, FCalcMethod.TAAM)
        w2_time = perf_counter() - start

        assert (
            w1_time / w2_time > 10
        ), "Atom assignment should be at least 10x slower than cache lookup"

        assert all(a == b for a, b in zip(w1.f_calc(2), w2.f_calc(2)))
        assert w1 is w2

    def test_different_structures(self, random_structure, tyrosine):
        DiscambWrapperCached.__cache = {}
        w1 = DiscambWrapperCached(random_structure)
        w2 = DiscambWrapperCached(tyrosine)
        assert w1 is not w2

    def test_different_unit_cell(self, random_structure, random_structure_u_iso):
        DiscambWrapperCached.__cache = {}

        assert "".join(s.label for s in random_structure.scatterers()) == "".join(
            s.label for s in random_structure_u_iso.scatterers()
        )
        assert random_structure_u_iso.unit_cell() != random_structure.unit_cell()

        w1 = DiscambWrapperCached(random_structure)
        w2 = DiscambWrapperCached(random_structure_u_iso)
        assert w1 is not w2

        random_structure.shake_sites_in_place(0.1)
        w2 = DiscambWrapperCached(random_structure)
        assert w1 is w2


class TestFromFile:
    def test_pdb(self, tmp_path, random_structure):
        pdb_file = tmp_path / "test.pdb"
        with pdb_file.open("w") as f:
            f.write(random_structure.as_pdb_file())

        DiscambWrapper.from_file(pdb_file)
        DiscambWrapper.from_file(str(pdb_file))

    def test_file_errors(self, tmp_path, random_structure):
        # Write a file with actual structure content
        pdb_file = tmp_path / "test.tmp"
        with pdb_file.open("w") as f:
            f.write(random_structure.as_pdb_file())

        with pytest.raises(
            ValueError, match="Supported files are .cif, .mmcif and .pdb. Got .tmp"
        ):
            DiscambWrapper.from_file(pdb_file)

        with pytest.raises(FileNotFoundError):
            DiscambWrapper.from_file(pdb_file.with_name("not_found"))

    def test_cif(self, tmp_path, random_structure):
        cif_file = tmp_path / "tmp.cif"
        with cif_file.open("w") as f:
            random_structure.as_cif_simple(out=f, format="corecif")

        DiscambWrapper.from_file(cif_file)

    def test_cif_multiple_structures(self, tmp_path, random_structure):
        import iotbx.cif

        cif = iotbx.cif.model.cif()
        cif["A"] = random_structure.as_cif_block(format="corecif")
        cif["B"] = random_structure.as_cif_block(format="corecif")
        cif_file = tmp_path / "tmp.cif"
        with cif_file.open("w") as f:
            f.write(str(cif))

        with pytest.raises(ValueError, match="Multiple structures found in file"):
            DiscambWrapper.from_file(cif_file)


class TestFromPdbCode:
    def test_working(self):
        DiscambWrapper.from_pdb_code("1HUF")

    def test_wrong_length(self):
        with pytest.raises(ValueError, match="pdb code must be 4 characters long"):
            DiscambWrapper.from_pdb_code("12345")
        with pytest.raises(ValueError, match="pdb code must be 4 characters long"):
            DiscambWrapper.from_pdb_code("123")

    def test_illegal_characters(self):
        with pytest.raises(ValueError, match="pdb code must be alphanumeric"):
            DiscambWrapper.from_pdb_code("123!")

    def test_not_found(self):
        # Warning: might fail if it exists in the future?
        with pytest.raises(ValueError, match="pdb code not found on rcsb.org"):
            DiscambWrapper.from_pdb_code("0000")
