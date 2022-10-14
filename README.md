# Ion-Transport-Model
The Ion Transport model accounts for the contributions of diffusion,convection and migration.This takes into account the Nernst-Planck equation, Fick's law of diffusion coupled with electroneutrality, and the Donnan equilibrium conditions.
In this study we use the finite element method (FEM) to find the numerical solutions of the differential equations developed from these conditions.

The simulations were preformed using the MATLAB bvp4c method.bvp4c is a finite difference code that implements the three-stage Lobatto IIIa formula.This is a collocation formula and the collocation polynomial provides a $C^1$-continuous solution that is fourth-order accurate uniformly in the interval of integration. Mesh selection and error control are based on the residual of the continuous solution.

The Concentration Equation:

![equation](https://latex.codecogs.com/svg.image?%5Cfrac%7Bd%5E%7B2%7D%7D%7Bd%20x%5E%7B2%7D%7DC_%7Bi%7D%5C!%20%5Cleft(x%5Cright)=-%5Cfrac%7Bz_%7Bi%7D%20%5Cleft(%5Cfrac%7Bd%7D%7Bd%20x%7DC_%7Bi%7D%5C!%20%5Cleft(x%5Cright)%5Cright)%20F%20%5Cleft(%5Cfrac%7Bd%7D%7Bd%20x%7D%5Cphi%20%5C!%20%5Cleft(x%5Cright)%5Cright)%7D%7BR%20T%7D-%5Cfrac%7Bz_%7Bi%7D%20C_%7Bi%7D%5C!%20%5Cleft(x%5Cright)%20F%20%5Cleft(%5Cfrac%7Bd%5E%7B2%7D%7D%7Bd%20x%5E%7B2%7D%7D%5Cphi%20%5C!%20%5Cleft(x%5Cright)%5Cright)%7D%7BR%20T%7D&plus;%5Cfrac%7B%5Cleft(%5Cfrac%7Bd%7D%7Bd%20x%7DC_%7Bi%7D%5C!%20%5Cleft(x%5Cright)%5Cright)%20%5Cnu%7D%7B%5Cmathrm%7BD%7D_%7Bi%7D%7D)

The Potential Equation:

![equation](https://latex.codecogs.com/svg.image?%5Cfrac%7Bd%7D%7Bd%20x%7D%5Cphi%20%5C!%20%5Cleft(x%5Cright)=-%5Cfrac%7BR%20T%20%5Cleft(z_%7B%5Cmathit%7Bpos%7D%7D%20%5Cmathrm%7BD%7D_%7B%5Cmathit%7Bpos%7D%7D%20%5Cleft(%5Cfrac%7Bd%7D%7Bd%20x%7DC_%7B%5Cmathit%7Bpos%7D%7D%5C!%20%5Cleft(x%5Cright)%5Cright)&plus;z_%7B%5Cmathit%7Bneg%7D%7D%20%5Cmathrm%7BD%7D_%7B%5Cmathit%7Bneg%7D%7D%20%5Cleft(%5Cfrac%7Bd%7D%7Bd%20x%7DC_%7B%5Cmathit%7Bneg%7D%7D%5C!%20%5Cleft(x%5Cright)%5Cright)&plus;%5Cfrac%7Bi%7D%7BF%7D-%5Cleft(z_%7B%5Cmathit%7Bneg%7D%7D%20C_%7B%5Cmathit%7Bneg%7D%7D&plus;z_%7B%5Cmathit%7Bpos%7D%7D%20C_%7B%5Cmathit%7Bpos%7D%7D%5Cright)%20%5Cnu%20%5Cright)%7D%7BF%20%5Cleft(z_%7B%5Cmathit%7Bneg%7D%7D%5E%7B2%7D%20%5Cmathrm%7BD%7D_%7B%5Cmathit%7Bneg%7D%7D%20C_%7B%5Cmathit%7Bneg%7D%7D&plus;z_%7B%5Cmathit%7Bpos%7D%7D%5E%7B2%7D%20%5Cmathrm%7BD%7D_%7B%5Cmathit%7Bpos%7D%7D%20C_%7B%5Cmathit%7Bpos%7D%7D%5Cright)%7D)
