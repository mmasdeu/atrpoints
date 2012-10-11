atrpoints
=========

Numerical verification of conjecture of Darmon on ATR points.

The main files to look at are tests.sage and tests_withlevel.sage. The files are commented and contain examples of the type of calculations that can be made. These examples use the class TestDL, which encapsulates the main classes of the code: HilbertModularForm and ATRField.

The files atrpoints.sage and atrpoints_withlevel.sage contain the code for the classes TestDL, HilbertModularForm, ATRField and for other auxiliary classes and methods. We used cython for the most computationally demanding functions (e.g., for computing Fourier coefficients of Hilbert modular forms and for computing their integrals): this is the content of atr_cython.spyx and atr_cython_withlevel.spyx.

The file quadratic_continued_fraction.sage contains the code for computing continued fractions of real quadratic fields. It is based on the sage patch http://trac.sagemath.org/sage_trac/ticket/11380, but with some modifications that allow for the computation of many different continued fractions of the elements in the field.
