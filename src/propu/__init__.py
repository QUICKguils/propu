"""Aerospace Propulsion computer code.

This python package holds all the computer code written
as part of the Aerospace Propulsion course (AERO0014).

Subpackages
-----------
exercices -- Collection of solved exercises.
project   -- Code developed as part of the project.
exam      -- Code developed as part of the written exam.

Submodules
----------
mplrc        -- Set some global matplotlib parameters.
constant     -- Constant quantities used throughout the code.
iteralg      -- Iterative algorithm utilities.
bemt         -- Blade element momentum techniques.
turbine      -- Gas turbine analyses.
"""

import pathlib

from propu import mplrc

_ROOT_PATH = pathlib.Path(__file__).parent.parent.parent
_SRC_PATH = pathlib.Path(__file__).parent

mplrc.load_rcparams()
