# Ion-Transport-Model
The Ion Transport model accounts for the contributions of diffusion,convection and migration.This takes into account the Nernst-Planck equation, Fick's law of diffusion coupled with electroneutrality, and the Donnan equilibrium conditions.
In this study we use the finite element method (FEM) to find the numerical solutions of the differential equations developed from these conditions.

The simulations were preformed using the MATLAB bvp4c method.bvp4c is a finite difference code that implements the three-stage Lobatto IIIa formula.This is a collocation formula and the collocation polynomial provides a $C^1$-continuous solution that is fourth-order accurate uniformly in the interval of integration. Mesh selection and error control are based on the residual of the continuous solution.
