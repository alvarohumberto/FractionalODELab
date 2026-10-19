(* ::Package:: *)

(* ::Title:: *)
(*FractionalODELab - Caputo fractional ODE solvers and Caputo-derivative operators*)


(* ::Text:: *)
(*Version 0.1.0. MIT License. Author: Alvaro H. Salas (ahsalass@unal.edu.co).*)
(*Problem: D_C^alpha y(t) = f(t, y(t)), y^(k)(0) = y0_k, k = 0..Ceiling[alpha]-1, 0 < alpha < 2, y scalar or vector.*)
(*Solvers: fractional Adams-Bashforth-Moulton PECE, implicit product-integration trapezoidal rule,*)
(*Lubich convolution quadrature (BDF1/BDF2 with starting weights), variable-order PECE.*)
(*Operators: L1 and Strong-FDM (A. H. Salas) approximations of the Caputo derivative; residual certificate.*)


BeginPackage["FractionalODELab`"];

FOLPECE::usage =
"FOLPECE[f, alpha, y0, {T, N}] solves D^alpha y = f[t, y] on [0, T] with N steps by the fractional \
Adams-Bashforth-Moulton predictor-corrector. Returns {t, y}. Options: \"CorrectorIterations\" (1), \"Implicit\" (False). \
For 1 < alpha < 2 give y0 = {y(0), y'(0)}.";
FOLTrapezoidal::usage =
"FOLTrapezoidal[f, alpha, y0, {T, N}] is the implicit product-integration trapezoidal rule (the PECE corrector solved by Newton).";
FOLLubich::usage =
"FOLLubich[f, alpha, y0, {T, N}] is the implicit Lubich convolution quadrature with BDF1 (\"Order\" -> 1) or BDF2 \
(\"Order\" -> 2, default) generating function and \"StartingWeights\" -> True.";
FOLLubichWeights::usage = "FOLLubichWeights[alpha, N, order] gives the convolution weights w_0..w_N of delta(z)^alpha.";
FOLVariableOrderPECE::usage =
"FOLVariableOrderPECE[f, alphaFun, y0, {T, N}] is the variable-order PECE with alpha = alphaFun[t] in (0, 1).";
FOLCaputoL1::usage = "FOLCaputoL1[u, alpha, h] gives the L1 approximations of D^alpha u at t_1..t_N from nodal values u_0..u_N (0<alpha<1).";
FOLCaputoStrongFDM::usage =
"FOLCaputoStrongFDM[u, alpha, h, N] gives the Strong-FDM (A. H. Salas) approximations of D^alpha u at t_1..t_N for a \
function u. Option \"InitialDerivative\" -> u^(p-1)(0), p = Ceiling[alpha] + 1 (Automatic: one-sided difference).";
FOLResidual::usage =
"FOLResidual[f, alpha, {t, y}] gives r_n = D^alpha y(t_n) - f[t_n, y_n], n = 1..N, for a scalar solution. \
Option \"Operator\" -> \"StrongFDM\" (default; quintic interpolant) or \"L1\".";
FOLCaputoPower::usage = "FOLCaputoPower[p, alpha, t] gives the exact Caputo derivative of t^p.";
FOLInterpolation::usage = "FOLInterpolation[{t, y}] gives an interpolating function of a computed solution.";
FOLVersion::usage = "FOLVersion gives the package version.";

FOLVariableOrderPECE::range = "alpha(t) = `1` at t = `2` is outside (0, 1).";


Begin["`Private`"];

FOLVersion = "0.1.0";

(* ---------- helpers ---------- *)
(* for alpha > 1 a list {y(0), y'(0)} of scalars is the usual input *)
icData[y0_, alpha_] := Module[{m = Max[1, Ceiling[alpha]]},
   Which[
    m == 1, {N[Flatten[{y0}]]},
    ListQ[y0] && Length[y0] == m, N[Flatten[{#}]] & /@ y0,
    True, PadRight[{N[Flatten[{y0}]]}, m, {0. N[Flatten[{y0}]]}]]];

taylorPart[t_List, y0_, alpha_] := Module[{ic = icData[y0, alpha]},
   Total[Table[Outer[Times, If[k == 0, ConstantArray[1., Length[t]], N[t]^k/k!], ic[[k + 1]]], {k, 0, Length[ic] - 1}]]];

wrapF[f_, d_] := If[d == 1, Function[{t, y}, N[Flatten[{f[t, First[y]]}]]], Function[{t, y}, N[Flatten[{f[t, y]}]]]];

newton[g_, x0_, tol_, maxit_] := Module[{x = N[x0], gx, J, dx, m = Length[x0], e},
   Do[
    gx = g[x];
    J = Transpose@Table[e = 10.^-7 (1 + Abs[x[[i]]]); (g[x + e UnitVector[m, i]] - gx)/e, {i, m}];
    dx = LinearSolve[J, -gx]; x = x + dx;
    If[Max[Abs[dx]] < tol (1 + Max[Abs[x]]), Break[]], {maxit}];
   x];

unwrap[t_, y_, d_] := {t, If[d == 1, y[[All, 1]], y]};

(* ---------- PECE / trapezoidal ---------- *)
Options[FOLPECE] = {"CorrectorIterations" -> 1, "Implicit" -> False, "Tolerance" -> 10.^-13, "MaxIterations" -> 50};
FOLPECE[f_, alpha_?NumericQ, y0_, {T_?NumericQ, nSteps_Integer?Positive}, OptionsPattern[]] :=
 Module[{a = N[alpha], h, t, Tp, d, F, y, fv, cP, cC, k, kA, kA1, b, aw, j, yp, hist, n, tol, mit},
  h = N[T]/nSteps; t = h Range[0, nSteps];
  Tp = taylorPart[t, y0, a]; d = Length[First[Tp]]; F = wrapF[f, d];
  tol = OptionValue["Tolerance"]; mit = OptionValue["MaxIterations"];
  y = ConstantArray[0., {nSteps + 1, d}]; fv = y;
  y[[1]] = Tp[[1]]; fv[[1]] = F[t[[1]], y[[1]]];
  cP = h^a/Gamma[a + 1]; cC = h^a/Gamma[a + 2];
  k = N[Range[0, nSteps + 1]]; kA = k^a; kA1 = k^(a + 1);
  Do[
   j = Range[0, n];
   b = kA[[n + 2 - j]] - kA[[n + 1 - j]];
   aw = kA1[[n - j + 3]] + kA1[[n - j + 1]] - 2 kA1[[n - j + 2]];
   aw[[1]] = n^(a + 1) - (n - a) (n + 1)^a;
   yp = Tp[[n + 2]] + cP (b . fv[[1 ;; n + 1]]);
   hist = Tp[[n + 2]] + cC (aw . fv[[1 ;; n + 1]]);
   If[TrueQ[OptionValue["Implicit"]],
    yp = newton[(# - hist - cC F[t[[n + 2]], #]) &, yp, tol, mit],
    Do[yp = hist + cC F[t[[n + 2]], yp], {OptionValue["CorrectorIterations"]}]];
   y[[n + 2]] = yp; fv[[n + 2]] = F[t[[n + 2]], yp],
   {n, 0, nSteps - 1}];
  unwrap[t, y, d]];

Options[FOLTrapezoidal] = {"Tolerance" -> 10.^-13, "MaxIterations" -> 50};
FOLTrapezoidal[f_, alpha_, y0_, {T_, nSteps_}, opts : OptionsPattern[]] :=
  FOLPECE[f, alpha, y0, {T, nSteps}, "Implicit" -> True, opts];

(* ---------- Lubich convolution quadrature ---------- *)
FOLLubichWeights[alpha_?NumericQ, nSteps_Integer, order_Integer : 2] :=
 Module[{a = N[alpha], g, q},
  g = FoldList[#1 (1 - (a + 1)/#2) &, 1., Range[nSteps]];
  If[order == 1, Return[g]];
  q = FoldList[#1 (1 - (a + 1)/#2)/3 &, 1., Range[nSteps]];
  1.5^a ListConvolve[g, q, {1, 1}, 0][[1 ;; nSteps + 1]]];

startingExponents[a_, order_, maxCount_ : 4] := Take[
   Union@Select[Flatten@Table[k + l a, {k, 0, order}, {l, 0, Floor[order/a] + 1}], 0 < # < order - 10^-12 &],
   UpTo[maxCount]];

Options[FOLLubich] = {"Order" -> 2, "StartingWeights" -> True, "Tolerance" -> 10.^-13, "MaxIterations" -> 50};
FOLLubich[f_, alpha_?NumericQ, y0_, {T_?NumericQ, nSteps_Integer?Positive}, OptionsPattern[]] :=
 Module[{a = N[alpha], ord = OptionValue["Order"], h, t, Tp, d, F, w, ha, qs, s, W, js, Vinv, kk, R, u, res, gs,
   tol = OptionValue["Tolerance"], mit = OptionValue["MaxIterations"], n},
  h = N[T]/nSteps; t = h Range[0, nSteps];
  Tp = taylorPart[t, y0, a]; d = Length[First[Tp]]; F = wrapF[f, d];
  w = FOLLubichWeights[a, nSteps, ord]; ha = h^a;
  qs = If[TrueQ[OptionValue["StartingWeights"]], startingExponents[a, ord], {}];
  s = Min[Length[qs], nSteps];
  W = ConstantArray[0., {nSteps + 1, Max[s, 1]}];
  If[s > 0,
   js = N[Range[s]];
   Vinv = Inverse[Table[js^q, {q, qs[[1 ;; s]]}]];
   kk = N[Range[0, nSteps]];
   Do[
    R = Table[Total[Join[{Gamma[q + 1]/Gamma[q + 1 - a] n^(q - a)}, -Reverse[w[[1 ;; n]]] kk[[2 ;; n + 1]]^q]],
      {q, qs[[1 ;; s]]}];
    W[[n + 1, 1 ;; s]] = Vinv . R, {n, 1, nSteps}]];
  u = ConstantArray[0., {nSteps + 1, d}];
  res[nn_, un_, uu_] := Module[{acc},
    acc = w[[1]] un + If[nn > 1, Reverse[w[[2 ;; nn]]] . uu[[2 ;; nn]], 0.];
    If[s > 0, acc = acc + W[[nn + 1, 1 ;; s]] . uu[[2 ;; s + 1]]];
    acc - ha F[t[[nn + 1]], Tp[[nn + 1]] + un]];
  If[s > 0,
   gs = Function[{x}, Module[{uu = u}, uu[[2 ;; s + 1]] = Partition[x, d];
      Flatten@Table[res[nn, uu[[nn + 1]], uu], {nn, 1, s}]]];
   u[[2 ;; s + 1]] = Partition[newton[gs, ConstantArray[0., s d], tol, mit], d]];
  Do[u[[n + 1]] = newton[res[n, #, u] &, u[[n]], tol, mit], {n, s + 1, nSteps}];
  unwrap[t, u + Tp, d]];

(* ---------- variable order (A. H. Salas) ---------- *)
Options[FOLVariableOrderPECE] = {"CorrectorIterations" -> 1};
FOLVariableOrderPECE[f_, alphaFun_, y0_, {T_?NumericQ, nSteps_Integer?Positive}, OptionsPattern[]] :=
 Module[{h, t, y0v, d, F, y, fv, al, j, b, aw, yp, hist, n},
  h = N[T]/nSteps; t = h Range[0, nSteps];
  y0v = N[Flatten[{y0}]]; d = Length[y0v]; F = wrapF[f, d];
  y = ConstantArray[0., {nSteps + 1, d}]; fv = y;
  y[[1]] = y0v; fv[[1]] = F[0., y0v];
  Do[
   al = N[alphaFun[t[[n + 2]]]];
   If[! (0 < al < 1), Message[FOLVariableOrderPECE::range, al, t[[n + 2]]]; Return[$Failed, Module]];
   j = N[Range[0, n]];
   b = (n + 1 - j)^al - (n - j)^al;
   aw = (n - j + 2)^(al + 1) + (n - j)^(al + 1) - 2 (n - j + 1)^(al + 1);
   aw[[1]] = n^(al + 1) - (n - al) (n + 1)^al;
   yp = y0v + h^al/Gamma[al + 1] (b . fv[[1 ;; n + 1]]);
   hist = y0v + h^al/Gamma[al + 2] (aw . fv[[1 ;; n + 1]]);
   Do[yp = hist + h^al/Gamma[al + 2] F[t[[n + 2]], yp], {OptionValue["CorrectorIterations"]}];
   y[[n + 2]] = yp; fv[[n + 2]] = F[t[[n + 2]], yp],
   {n, 0, nSteps - 1}];
  unwrap[t, y, d]];

(* ---------- Caputo operators ---------- *)
FOLCaputoPower[p_, alpha_, t_] := Gamma[p + 1]/Gamma[p + 1 - alpha] t^(p - alpha);

FOLCaputoL1[u_List, alpha_?NumericQ, h_?NumericQ] := Module[{a = N[alpha], nn = Length[u] - 1, k, b, du, c},
   k = N[Range[0, nn]]; b = (k + 1)^(1 - a) - k^(1 - a); du = Differences[N[u]];
   c = h^-a/Gamma[2 - a];
   Table[c Reverse[b[[1 ;; n]]] . du[[1 ;; n]], {n, 1, nn}]];

Options[FOLCaputoStrongFDM] = {"InitialDerivative" -> Automatic};
FOLCaputoStrongFDM[u_, alpha_?NumericQ, h_?NumericQ, nSteps_Integer, OptionsPattern[]] :=
 Module[{a = N[alpha], p, d0, e = 10.^-4, jj, D, c, wgt},
  p = Ceiling[a] + 1;
  d0 = Replace[OptionValue["InitialDerivative"],
    Automatic :> If[p == 2, (u[e] - u[0.])/e, (u[2 e] - 2 u[e] + u[0.])/e^2]];
  jj = Range[nSteps];
  D = Sum[(-1)^(p - r) Binomial[p, r] (u[(# - 1 + r/p) h] & /@ jj), {r, 0, p}];
  c = p^p/(h^p Gamma[p + 1 - a]);
  Table[
   wgt = ((n - jj[[1 ;; n]] + 1) h)^(p - a) - ((n - jj[[1 ;; n]]) h)^(p - a);
   (n h)^(p - 1 - a)/Gamma[p - a] d0 + c wgt . D[[1 ;; n]], {n, 1, nSteps}]];

Options[FOLResidual] = {"Operator" -> "StrongFDM", "InitialDerivative" -> Automatic};
FOLResidual[f_, alpha_?NumericQ, {t_List, y_List}, OptionsPattern[]] :=
 Module[{h = t[[2]] - t[[1]], nn = Length[t] - 1, D, ip},
  D = If[OptionValue["Operator"] === "L1", FOLCaputoL1[y, alpha, h],
    ip = Interpolation[Transpose[{t, y}], InterpolationOrder -> 5];
    FOLCaputoStrongFDM[ip, alpha, h, nn, "InitialDerivative" -> OptionValue["InitialDerivative"]]];
  D - MapThread[N[f[#1, #2]] &, {Rest[t], Rest[y]}]];

FOLInterpolation[{t_List, y_List}] := If[MatrixQ[y],
   Interpolation[Transpose[{t, #}]] & /@ Transpose[y], Interpolation[Transpose[{t, y}]]];

End[];
EndPackage[];
