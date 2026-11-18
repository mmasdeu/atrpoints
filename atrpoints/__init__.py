r"""
atrpoints: ATR points over real quadratic fields.

Modules:

- :mod:`atrpoints.arithgrouprq` -- the groups SL_2(O_F) and Gamma_1(N) for a real
  quadratic field F, with presentations and word problem (:class:`ArithGroupRQ`).
- :mod:`atrpoints.quadcontfrac` -- continued fractions in real quadratic fields.
- :mod:`atrpoints.atrpoints` -- ATR point computations at level 1 (:class:`TestDL`).
- :mod:`atrpoints.atrpoints_withlevel` -- ATR point computations with level.
- :mod:`atrpoints.utility` -- auxiliary functions for the computations with level.

The two ATR modules define classes with the same names (``TestDL``, ``ATRField``, ...),
so import the one you need, e.g. ``from atrpoints.atrpoints_withlevel import *``.
"""
from .quadcontfrac import QuadraticContinuedFraction, all_qcf_up_to_length, quadratic_continued_fraction
from .arithgrouprq import ArithGroupRQ
