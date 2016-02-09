FE Solver
=========

* Author: Alex Griessman (<alex.griessman@gmail.com>)
* Repository <https://github.com/mirrorscotty/fe-solver/>

This is a basic framework for solving differential equations using finite
element method. Currently it is capable of solving both 1D and 2D steady-state
problems and transient problems.

Contents
--------
* `output` - Functions for outputting simulation data
* `scaling` - Nondemensionalization stuff for heat and mass transfer problems
* `solver` - Actual finite element solver (both linear and nonlinear)
  * `integration` - Three and five point Gaussian quadrature
  * `mesh` - Basic 1D and 2D mesh routines
  * `ode` - Simple ODE solver using Euler's method

Building
--------
This code was designed to be compiled as a library and statically linked. To build just the library, execute `make fe-solver.a`. Documentation can be built using Doxygen using `make doc`.

Usage
-----
Include `fe-solver.h` in any C file that requires functions provided by this library, and instruct the compiler to link the executable against `fe-solver.a`. The PDE's to be solved should be written in weak form and then differentiated with respect to each dependent variable. Three functions should be written for each derivative, the time-dependent portion, the time-independent portion, and the load vector. Examples are provided [here](https://github.com/mirrorscotty/fe-problems).

Dependencies
------------
Compilation requires GCC and GNU Make. Building documentation requires Doxygen and LaTeX to generate PDF output. This code also requires the matrix library found [here](https://github.com/mirrorscotty/matrix).
