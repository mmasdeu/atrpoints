atrpoints
=========

Numerical verification of Darmon's conjecture on ATR points, for elliptic curves over
real quadratic fields, together with the arithmetic groups SL_2(O_F) and Gamma_1(N)
(presentations and word problem) needed for the computations with level.

Requires [SageMath](https://www.sagemath.org) and, for `atrpoints.arithgrouprq`,
[darmonpoints](https://github.com/mmasdeu/darmonpoints).

Installation
------------

From the repository directory, in the Python environment of your Sage installation:

    sage -pip install --no-build-isolation -e .

This compiles the Cython extensions (`atrpoints/atr_cython*.pyx`). To use the code without
installing it, build them in place with `sage -python setup.py build_ext --inplace` and
add the repository directory to `PYTHONPATH`.

Package contents
----------------

- `atrpoints/arithgrouprq.py` -- `ArithGroupRQ(F, level=N)`: the group Gamma_1(N) in
  SL_2(O_F) (SL_2(O_F) itself for N = 1), with a presentation computed by
  Reidemeister-Schreier and a word-problem solver. It is a darmonpoints arithmetic group,
  so abelianizations, bounding cycles, etc. work on it. See `IMPLEMENTATION.md`.
- `atrpoints/quadcontfrac.py` -- continued fractions in real quadratic fields, based on
  the Sage patch http://trac.sagemath.org/sage_trac/ticket/11380, modified to compute
  many different continued fractions of the same element.
- `atrpoints/atrpoints.py` -- ATR points at level 1: the classes `TestDL`,
  `HilbertModularForm`, `ATRField` and auxiliary functions.
- `atrpoints/atrpoints_withlevel.py` and `atrpoints/utility.py` -- the same for curves with
  nontrivial conductor.
- `atrpoints/atr_cython.pyx`, `atrpoints/atr_cython_withlevel.pyx` -- Cython code for the
  most demanding parts: Fourier coefficients of the Hilbert modular forms and the integrals.

The two ATR modules define classes with the same names, so import the one you need:

    from atrpoints.atrpoints import *            # level 1
    from atrpoints.atrpoints_withlevel import *  # with level

Examples
--------

The scripts in `examples/` are commented and show the kind of computations that can be made:

- `atr_example.sage` -- checking the conjecture with `TestDL` (level 1);
- `main.sage` -- computations with `ArithGroupRQ`.

Run them with `sage examples/<name>.sage`.

Tests
-----

    sage -python -m pytest tests
