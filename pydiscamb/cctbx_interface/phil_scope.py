from typing import Any

import iotbx.phil
from libtbx.phil import scope_extract

from pydiscamb.taam_parameters import get_default_databank

pydiscamb_master_params = iotbx.phil.parse(
    f"""
    taam 
    .help = Transferrable Aspherical Atom Model (TAAM)-specific parameters
    {{
        unit_cell_charge = 0.0
            .type = float
            .help = Unit cell charge, used to scale multipolar parameters
        scale_pval_to_charge = True
            .type = bool
            .help = Whether to scale multipolar parameters to fit the given charge
        nproc = 1
            .type = int
            .help = Number of cores to use for computing
        freeze_local_coordinate_system = False
            .type = bool
            .help = Whether to re-calculate local coordinate systems for all atoms when updating the structure
        bank_path = {get_default_databank()}
            .type = str
            .help = Path to TAAM parameter bank file
            .expert_level = 2
    }}
"""
)


def scope_to_taam_dict(scope: scope_extract) -> dict[str, Any]:
    assert hasattr(scope, "discamb")
    taam_params = scope.discamb.taam
    out = {
        key: getattr(taam_params, key)
        for key in dir(taam_params)
        if key[:2] != "__" and getattr(taam_params, key) is not None
    }
    # Manually translate params renamed in the scope
    out["scale"] = out.pop("scale_pval_to_charge")
    out["n_cores"] = out.pop("nproc")
    out["frozen_lcs"] = out.pop("freeze_local_coordinate_system")

    return out
