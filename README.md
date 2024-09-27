# TAAM_SF

Simple pybind11 wrapper to communicate with discamb from cctbx

## Installation

Only tested with Cmake 3.22.3, gcc 11.4 on WSL2 ubuntu in Windows 11.
These can also be installed with conda (`conda install gcc=11.4`).

```bash
conda create --name taam_sf_dev python=3.8 -y
conda activate taam_sf_dev
conda install openmp -c conda-forge -y
pip install .
```

## Build with conda
Install conda-build [(recommended to do in base env)](https://docs.conda.io/projects/conda-build/en/latest/install-conda-build.html)

```bash
conda activate base
conda install conda-build
```

Build the package
```bash
cd TAAM_SF
conda build ./conda -c conda-forge
```

## Example usage

```python
from __future__ import absolute_import, division, print_function

import taam_sf

def run():
  # Generate a random structure
  from cctbx.development import random_structure as cctbx_random_structure
  from cctbx.sgtbx import space_group_info

  group = space_group_info(19)
  xrs = cctbx_random_structure.xray_structure(
      space_group_info=group,
      elements=["Ni", "C"] * 3,
  )
  xrs.scattering_type_registry(table = "electron")

  d_min = 2

  # Calculate structure factors with cctbx
  fcalc= xrs.structure_factors(d_min=d_min).f_calc()
  fcalc_cctbx = fcalc.data()

  # Calculate structure factors with DiSCaMB
  fcalc_discamb = taam_sf.calculate_structure_factors_IAM(xrs, d_min)

  # Output the difference of the results
  diff = [fc - fd for fc, fd in zip(fcalc_cctbx, fcalc_discamb)]
  print(f"Number of reflections: {len(diff)}")
  print(f"Mean f_calc deviation: {sum(diff) / len(diff) :.1e}")

if (__name__ == "__main__"):
  run()
```

## Testing

```bash
conda install cctbx-base
pip install pytest
pytest
```
