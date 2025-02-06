"""
propu -- Aerospace Propulsion computer code
===========================================

This python package holds all the computer code written
as part of the Aerospace Propulsion course (AERO0014).

Subpackages
-----------
util      -- Common utilities that are used throughout the code.
bemt      -- Blade element momentum theory and associated project.
turbine   -- Gas turbine analyses and associated project.
exercices -- Collection of solved exercises.
exam      -- Code developed as part of the written exam.
"""

import pathlib as _pathlib

from propu.util import mplrc

ROOT_PATH = _pathlib.Path(__file__).parent.parent.parent

mplrc.load_rcparams()
