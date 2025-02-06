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
source .venv/bin/activate  # on Unix/macOS
venv\Scripts\activate      # on Windows
```

Still from the top level, install the package.
```sh
python -m pip install .
```

## Usage

```python
import propu
# ...do whathever you want with the package
```

## Project architecture

- `src/`
  - `util\` Common utilities that are used throughout the code.
    - `mplrc.py` Set some global matplotlib parameters.
    - `constant.py` Constant quantities used throughout the code.
    - `isatmosphere.py` International Standard Atmosphere.
  - `bemt\` Blade element momentum theory and associated project.
  - `turbine\` Gas turbine analyses and associated project.
  - `exercices\` Collection of solved exercises.
  - `exam\` Code developed as part of the written exam.
