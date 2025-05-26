# Aerospace propulsion

This python package holds all the computer code written
as part of the Aerospace Propulsion course (AERO0014).

## Installation

From the top level directory (where this README lies),
create a virtual environment, and activate it.
```sh
# create a virtual environment
python -m venv venv

# activate the virtual environment
source venv/bin/activate  # on Unix/macOS
venv\Scripts\activate     # on Windows
```

Still from the top level, install the package
(optionally in editable mode, with the `-e` flag).
```sh
python -m pip install -e .
```

## Usage

```python
import propu
# ...do whathever you want with the package

# Example 1
# Run the second part of the project.
from propu.project import part2
sol = part2.main()

# Example 2
# Import all solved execises (recommended in interactive session only),
# and execute one of them.
from propu.exercises import *
E4_2__Nozzle_change.main()

# Example 3
# Get the air density at 30_000 ft.
from propu.constant import mconv, get_isa
rho = get_isa(mconv(30_000, "ft", "m")).rho
```

## Project layout

The source code lies in `src/propu/`.
It contains the following packages and modules.
- Packages:
  - `exercises/` Collection of solved exercises.
  - `project/` Code developed as part of the project.
  - `exam/` Code developed as part of the written exam.
- Modules:
  - `mplrc.py` Set some global Matplotlib parameters.
  - `constant.py` Constant quantities used throughout the code.
  - `iteralg.py` Iterative algorithm utilities.
  - `bemt.py` Blade element momentum techniques.
  - `turbine.py` Gas turbine and jet engine analyses.
