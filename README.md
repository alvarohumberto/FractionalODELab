# FractionalODELab

Solvers for Caputo fractional differential equations  D^alpha y = f(t, y),  0 < alpha < 2
(scalar equations and commensurate systems), in the Wolfram Language, with a cross-validated Python port.

* `FOLPECE`: fractional Adams-Bashforth-Moulton predictor-corrector (Diethelm-Ford-Freed), any number of corrector iterations
* `FOLTrapezoidal`: implicit product-integration trapezoidal rule (Newton)
* `FOLLubich`: Lubich convolution quadrature, BDF1/BDF2, with starting weights
* `FOLVariableOrderPECE`: variable-order PECE for alpha = alpha(t)
* `FOLCaputoL1`, `FOLCaputoStrongFDM`: discrete Caputo operators
* `FOLResidual`: a-posteriori residual certificate of a computed solution

## Quick start (Mathematica / Wolfram Engine >= 14, tested with 15.0)
```wolfram
Get["mathematica/FractionalODELab.wl"]
{t, y} = FOLPECE[-#2 &, 0.6, 1., {5., 500}, "CorrectorIterations" -> 2];
Max[Abs[y - MittagLefflerE[0.6, -t^0.6]]]
```
For 1 < alpha < 2 give `y0 = {y(0), y'(0)}`; for systems give vectors.

## Reproduce the paper
```
wolframscript -file mathematica/RunTests.wls      # 36 tests
wolframscript -file mathematica/MakeFigures.wls   # figures/ and tables/
```
Python (optional): `pip install -e .` then `pytest tests` (20 tests).

## Layout
`mathematica/` package, tests, figure script · `src/fractional_ode_lab/` Python port · `tests/` Python tests ·
`benchmarks/` Python reference values used for cross-validation · `figures/`, `tables/` outputs · `paper/` manuscript.

## License
MIT. Author: Alvaro H. Salas, Universidad Nacional de Colombia (ahsalass@unal.edu.co).
