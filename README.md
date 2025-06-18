# pyDiSCaMB

![Tests](https://github.com/viljarjf/pyDiSCaMB/actions/workflows/test.yaml/badge.svg?branch=main)

Simple pybind11 wrapper to communicate with DiSCaMB from cctbx

`git clone --recursive git@github.com:viljarjf/pyDiSCaMB.git`

**LICENSE NOTE**
Recursive cloning installs the [MATTS databank](https://www.github.com/discamb-project/MATTS), which uses a license restricting commercial use.
Building the project with the databank will **include the databank in the built module**, which should be taken into account if the build is distributed.

## Installation

Tested with Cmake 3.22.3, gcc 11.4 on WSL2 ubuntu in Windows 11.
These can be installed with conda (`conda install gcc=11.4`).
Also verified using MSVC 19.41.

```bash
conda create --name pyDiSCaMB_dev python=3.8 -y
conda activate pyDiSCaMB_dev
conda install openmp -c conda-forge -y
pip install .
```

## Build with conda
Tested on Ubuntu WSL2.

Install conda-build [(recommended to do in base env)](https://docs.conda.io/projects/conda-build/en/latest/install-conda-build.html)

```bash
conda activate base
conda install conda-build
```

Build the package
```bash
cd pyDiSCaMB
conda build ./conda -c conda-forge
```

## Example usage

```python
from __future__ import absolute_import, division, print_function

import pydiscamb

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
  fcalc_discamb = pydiscamb.calculate_structure_factors_IAM(xrs, d_min)

  # Output the difference of the results
  diff = [fc - fd for fc, fd in zip(fcalc_cctbx, fcalc_discamb)]
  print(f"Number of reflections: {len(diff)}")
  print(f"Mean f_calc deviation: {sum(diff) / len(diff) :.1e}")

if (__name__ == "__main__"):
  run()
```

## Testing

```bash
conda install cctbx-base pytest -c conda-forge
pytest
```
