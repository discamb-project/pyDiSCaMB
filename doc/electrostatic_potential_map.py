import pydiscamb
from cctbx import miller
from iotbx import pdb
from pathlib import Path

# Replace this with a path to your structure file
structure_file = Path(__file__).resolve().parent / "data" / "tyrosine.pdb"

# Define resolution
d_min = 2.5

# Load structure
xrs = pdb.input(file_name=str(structure_file)).xray_structure_simple()
xrs.scattering_type_registry(table="electron")

# Prepare calculators
taam_calculator = pydiscamb.DiscambWrapper(
    xrs,
    method=pydiscamb.FCalcMethod.TAAM,
)

# Log the number of typed atoms
log = open("epm_atom_typing.log", "w")
taam_calculator.show_atom_type_assignment(log)
log.close()

# Set resolution
miller_set = miller.build_set(
    crystal_symmetry=xrs.crystal_symmetry(),
    anomalous_flag=False,
    d_min=d_min,
)
taam_calculator.set_indices(miller_set.indices())

# Compute structure factors
fcalc_taam = taam_calculator.f_calc()

# Create a miller.array to save
electrostatic_potential_map_array = miller_set.array(
    data=fcalc_taam,
)
electrostatic_potential_map_array.write_mtz(
    file_name="electrostatic_potential_map.mtz",
)
