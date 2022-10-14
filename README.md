# Ion-Transport-Model
The Ion Transport model accounts for the contributions of diffusion,convection and migration.This takes into account the Nernst-Planck equation, Fick's law of diffusion coupled with electroneutrality, and the Donnan equilibrium conditions.
In this study we use the finite element method (FEM) to find the numerical solutions of the differential equations developed from these conditions.

The simulations were preformed using the MATLAB bvp4c method.bvp4c is a finite difference code that implements the three-stage Lobatto IIIa formula.This is a collocation formula and the collocation polynomial provides a $C^1$-continuous solution that is fourth-order accurate uniformly in the interval of integration. Mesh selection and error control are based on the residual of the continuous solution.

The Concentration Equation:

$\frac{d^{2}}{d x^{2}}C_{i} \left(x\right)=-\frac{z_{i} \left(\frac{d}{d x}C_{i} \left(x\right)\right) F \left(\frac{d}{d x}\phi  \left(x\right)\right)}{R T}-\frac{z_{i} C_{i} \left(x\right) F \left(\frac{d^{2}}{d x^{2}}\phi \! \left(x\right)\right)}{R T}+\frac{\left(\frac{d}{d x}C_{i} \left(x\right)\right) \nu}{\mathrm{D}_{i}}$

The Potential Equation:


https://latex.codecogs.com/svg.image?\frac{d}{d&space;x}\phi&space;\!&space;\left(x\right)=-\frac{R&space;T&space;\left(z_{\mathit{pos}}&space;\mathrm{D}_{\mathit{pos}}&space;\left(\frac{d}{d&space;x}C_{\mathit{pos}}\!&space;\left(x\right)\right)&plus;z_{\mathit{neg}}&space;\mathrm{D}_{\mathit{neg}}&space;\left(\frac{d}{d&space;x}C_{\mathit{neg}}\!&space;\left(x\right)\right)&plus;\frac{i}{F}-\left(z_{\mathit{neg}}&space;C_{\mathit{neg}}&plus;z_{\mathit{pos}}&space;C_{\mathit{pos}}\right)&space;\nu&space;\right)}{F&space;\left(z_{\mathit{neg}}^{2}&space;\mathrm{D}_{\mathit{neg}}&space;C_{\mathit{neg}}&plus;z_{\mathit{pos}}^{2}&space;\mathrm{D}_{\mathit{pos}}&space;C_{\mathit{pos}}\right)}
