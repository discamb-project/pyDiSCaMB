## Test script that will open 1000 cached wrapper objects
# Import modules
from pydiscamb.discamb_wrapper import DiscambWrapperCached, FCalcMethod
from cctbx import crystal, xray
from cctbx.array_family import flex

# Define xrs
xrs = xray.structure(
    crystal_symmetry=crystal.symmetry(
        unit_cell=(10, 10, 10, 90, 90, 90),
        space_group_symbol="P1",
    ),
    scatterers=flex.xray_scatterer([xray.scatterer("C1", site=(0, 0, 0))]),
)
# Get a bunch of wrapper objects
for _ in range(1000):
    DiscambWrapperCached(xrs, FCalcMethod.TAAM)
