#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Verification suite for FractionalODELab (Wolfram Language).  Run: wolframscript -file RunTests.wls
   Writes test_report.txt / test_report.csv next to this script. *)

$dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{$dir, "FractionalODELab.wl"}]];
$root = ParentDirectory[$dir];

rate[{e1_, e2_}] := Log2[e1/e2];
ml[a_, t_] := N[MittagLefflerE[a, -t^a], 30];
bench[a_] := Function[{t, y}, 40320/Gamma[9 - a] t^(8 - a) - 3 Gamma[5 + a/2]/Gamma[5 - a/2] t^(4 - a/2) +
     9/4 Gamma[a + 1] + (3/2 t^(a/2) - t^4)^3 - Abs[y]^(3/2)];
benchExact[a_, t_] := t^8 - 3 t^(4 + a/2) + 9/4 t^a;
y0R[a_] := If[a > 1, {1., 0.}, 1.];
y0B[a_] := If[a > 1, {0., 0.}, 0.];
lastErr[sol_, ex_] := Abs[Last[sol[[2]]] - ex];

tests = Join[
  (* PECE on the Mittag-Leffler relaxation, observed order min(2, 1 + alpha) *)
  Table[With[{a = a}, VerificationTest[
     Module[{e = Table[lastErr[FOLPECE[-#2 &, a, y0R[a], {2., n}], ml[a, 2.]], {n, {200, 400}}]},
      {Last[e] < 2*^-5, rate[e] > Min[2, 1 + a] - 0.25}], {True, True}, TestID -> "PECE_MittagLeffler_a=" <> ToString[a]]],
   {a, {0.3, 0.5, 0.8, 1.0, 1.5, 1.8}}],
  (* second order on the Diethelm-Ford-Freed benchmark *)
  Flatten@Table[With[{a = a, m = m}, VerificationTest[
      Module[{e = Table[lastErr[
            Switch[m, "PECE3", FOLPECE[bench[a], a, y0B[a], {1., n}, "CorrectorIterations" -> 3],
             "TRAP", FOLTrapezoidal[bench[a], a, y0B[a], {1., n}],
             "LUB2", FOLLubich[bench[a], a, y0B[a], {1., n}]], benchExact[a, 1.]], {n, {160, 320}}]},
       {Last[e] < 2*^-4, rate[e] > 1.6}], {True, True}, TestID -> "Benchmark_" <> m <> "_a=" <> ToString[a]]],
    {a, {0.25, 0.5, 0.75, 1.25, 1.5, 1.75}}, {m, {"PECE3", "TRAP", "LUB2"}}],
  {
   (* system: fractional Lotka-Volterra *)
   VerificationTest[Module[{fLV = {#2[[1]] - 0.5 #2[[1]] #2[[2]], -#2[[2]] + 0.5 #2[[1]] #2[[2]]} &, yf, y},
      yf = Last[FOLPECE[fLV, 0.8, {2., 1.}, {5., 3200}][[2]]]; y = Last[FOLPECE[fLV, 0.8, {2., 1.}, {5., 400}][[2]]];
      Max[Abs[y - yf]] < 10^-3], True, TestID -> "System_LotkaVolterra"],
   (* variable order with constant alpha reproduces PECE exactly *)
   VerificationTest[Max[Abs[FOLVariableOrderPECE[-#2 &, 0.7 &, 1., {2., 200}][[2]] -
        FOLPECE[-#2 &, 0.7, 1., {2., 200}][[2]]]] < 10^-14, True, TestID -> "VariableOrder_constant"],
   (* variable order self-convergence *)
   VerificationTest[Module[{al = 0.75 + 0.2/(1 + Exp[-1.5 (# - 3)]) &, fl = 2 #2 (1 - #2) &, yr, e},
      yr = Last[FOLVariableOrderPECE[fl, al, 0.2, {10., 4000}][[2]]];
      e = Table[Abs[Last[FOLVariableOrderPECE[fl, al, 0.2, {10., n}][[2]]] - yr], {n, {250, 500}}];
      rate[e] > 1.7], True, TestID -> "VariableOrder_selfconvergence"],
   (* Lubich weights: BDF1 weights are the Grunwald-Letnikov coefficients *)
   VerificationTest[Max[Abs[FOLLubichWeights[0.6, 10, 1] - Table[(-1)^j Binomial[0.6, j], {j, 0, 10}]]] < 10^-14,
    True, TestID -> "LubichWeights_GL"]
   },
  (* Strong-FDM: second order for smooth functions, 0 < alpha < 2 *)
  Table[With[{a = a}, VerificationTest[
     Module[{ex = N[Sum[1/Gamma[k + 1 - a], {k, Ceiling[a], 80}], 20], e},
      e = Table[Abs[Last[FOLCaputoStrongFDM[Exp, a, 1./n, n, "InitialDerivative" -> 1.]] - ex], {n, {80, 160}}];
      rate[e] > 1.9], True, TestID -> "StrongFDM_exp_a=" <> ToString[a]]], {a, {0.3, 0.5, 0.8, 1.3, 1.7}}],
  {
   VerificationTest[Module[{a = 0.5, e}, e = Table[Abs[Last[FOLCaputoL1[N[Range[0, n]/n]^3.5, a, 1./n]] -
          FOLCaputoPower[3.5, a, 1.]], {n, {80, 160}}]; Abs[rate[e] - (2 - a)] < 0.1], True, TestID -> "L1_order"],
   VerificationTest[Module[{a = 0.5, f, r},
      f = Function[{t, y}, Gamma[4]/Gamma[4 - a] t^(3 - a) + 1 + t^3 - y];
      r = Table[Max[Abs[FOLResidual[f, a, FOLPECE[f, a, 1., {1., n}, "CorrectorIterations" -> 3], "InitialDerivative" -> 0.]]],
        {n, {80, 160}}]; r[[2]] < r[[1]]/3], True, TestID -> "Residual_decreases"],
   (* cross-implementation agreement with the Python package *)
   VerificationTest[Module[{rows = Rest@Import[FileNameJoin[{$root, "benchmarks", "python_reference.csv"}], "CSV"], v},
      Max@Table[
        v = Last[Switch[r[[2]],
            "PECE", FOLPECE[#1, r[[3]], #2, {#3, r[[4]]}],
            "PECE3", FOLPECE[#1, r[[3]], #2, {#3, r[[4]]}, "CorrectorIterations" -> 3],
            "TRAP", FOLTrapezoidal[#1, r[[3]], #2, {#3, r[[4]]}],
            "LUB2", FOLLubich[#1, r[[3]], #2, {#3, r[[4]]}]] &[
           If[r[[1]] == "relax", -#2 &, bench[r[[3]]]], If[r[[1]] == "relax", y0R[r[[3]]], y0B[r[[3]]]],
           If[r[[1]] == "relax", 2., 1.]][[2]]];
        Abs[v - r[[5]]]/(1 + Abs[r[[5]]]), {r, rows}] < 10^-10], True, TestID -> "Python_cross_validation"]
   }];

report = TestReport[tests];
lines = Table[{t["TestID"], ToString[t["Outcome"]], Round[1000 QuantityMagnitude[t["AbsoluteTimeUsed"], "Seconds"], 0.1]},
   {t, Values[report["TestResults"]]}];
Export[FileNameJoin[{$dir, "test_report.csv"}], Prepend[lines, {"test", "outcome", "ms"}], "CSV"];
summary = StringJoin["FractionalODELab (Wolfram Language) ", FOLVersion, " - ", $Version, "\n",
   "Passed ", ToString[report["TestsSucceededCount"]], " / ", ToString[Length[tests]], "\n\n",
   StringRiffle[StringRiffle[ToString /@ #, "  "] & /@ lines, "\n"], "\n"];
Export[FileNameJoin[{$dir, "test_report.txt"}], summary, "Text"];
Print[summary];
Do[If[t["Outcome"] =!= "Success", Print[t["TestID"], ": ", t["ActualOutput"], "  ", t["ActualMessages"]]],
  {t, Values[report["TestResults"]]}];
