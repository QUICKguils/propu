"""
util -- Common utilities that are used throughout the code
==========================================================

Submodules
----------
mplrc        -- Set some global matplotlib parameters.
constant     -- Constant quantities used throughout the code.
isatmosphere -- International Standard Atmosphere.

Note
----
Using ``from util import *`` will bring all the submodules into
the current namespace, under short aliases:
  mplrc        -> mplrc
  constant     -> c
  isatmosphere -> isa
"""

from . import mplrc
from . import constant as c
from . import isatmosphere as isa

print(dir())

__all__ = [s for s in dir() if not s.startswith('_')]
