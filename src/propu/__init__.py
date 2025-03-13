"""
Aerospace Propulsion computer code
==================================

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
isatmosphere -- International Standard Atmosphere.
bemt         -- Blade element momentum techniques.
turbine      -- Gas turbine analyses.
"""

import pathlib as _pathlib

from propu import mplrc

ROOT_PATH = _pathlib.Path(__file__).parent.parent.parent

mplrc.load_rcparams()
