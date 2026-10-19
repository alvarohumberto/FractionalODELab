#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Figures, tables and manuscript numbers for FractionalODELab.  Run: wolframscript -file MakeFigures.wls
   Output: ../figures/*.pdf|png, ../tables/*.tex|csv *)

$dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{$dir, "FractionalODELab.wl"}]];
$root = ParentDirectory[$dir];
$fig = FileNameJoin[{$root, "figures"}]; $tab = FileNameJoin[{$root, "tables"}];
Quiet[CreateDirectory /@ {$fig, $tab}];
$t0 = AbsoluteTime[];
log[s__] := Print[Round[AbsoluteTime[] - $t0, 0.1], "s  ", s];

(* ---------------- style ---------------- *)
cols = RGBColor /@ {"#1f5fa8", "#d1661a", "#2a8c55", "#7a4fa3", "#c0392b", "#5d6d7e", "#b8860b", "#17becf"};
font = {FontFamily -> "Times New Roman", FontSize -> 11};
frameOpts = Sequence[Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.7]],
   FrameTicksStyle -> Directive[Black, 10], BaseStyle -> font, ImageSize -> 330, AspectRatio -> 0.72,
   ImagePadding -> {{55, 12}, {40, 10}}];
lbl[g_, s_String] := Labeled[g, Style[s, Bold, 12, FontFamily -> "Times New Roman"], {{Top, Left}}];
row[gs__] := Grid[{{gs}}, Spacings -> {1.5, 0}, Alignment -> Top];
save[name_, g_] := (Export[FileNameJoin[{$fig, name <> ".pdf"}], g];
   Export[FileNameJoin[{$fig, name <> ".png"}], g, ImageResolution -> 300]; log["wrote ", name]);
leg[cs_, ls_, ___] := Placed[LineLegend[cs, ls, LegendMarkerSize -> {20, 10}, LabelStyle -> Append[font, 10],
    LegendLayout -> {"Row", 2}], Below];
mleg[cs_, ms_, ls_] := Placed[LineLegend[cs, ls, LegendMarkers -> ms, LegendMarkerSize -> {22, 10},
    LabelStyle -> Append[font, 10], LegendLayout -> {"Row", 2}], Below];
mk = {"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]", "\[FilledUpTriangle]", "\[EmptyCircle]", "\[EmptySquare]"};
slope[{x1_, x2_}, y1_, p_] := {{x1, y1}, {x2, y1 (x2/x1)^-p}};

bench[a_] := Function[{t, y}, 40320/Gamma[9 - a] t^(8 - a) - 3 Gamma[5 + a/2]/Gamma[5 - a/2] t^(4 - a/2) +
     9/4 Gamma[a + 1] + (3/2 t^(a/2) - t^4)^3 - Abs[y]^(3/2)];
benchExact[a_, t_] := t^8 - 3 t^(4 + a/2) + 9/4 t^a;
y0B[a_] := If[a > 1, {0., 0.}, 0.];

(* ---------------- Fig. 1: Mittag-Leffler relaxation ---------------- *)
relaxPanel[as_, cs_, T_, n_] := Module[{sols},
   sols = Table[FOLPECE[-#2 &, a, If[a > 1, {1., 0.}, 1.], {T, n}, "CorrectorIterations" -> 2], {a, as}];
   Show[
    Plot[Evaluate[Table[MittagLefflerE[a, -t^a], {a, as}]], {t, 0, T}, PlotStyle -> (Directive[#, AbsoluteThickness[1.4]] & /@ cs),
     PlotRange -> All],
    ListPlot[Table[Transpose[sols[[i]]][[1 ;; -1 ;; n/20]], {i, Length[as]}],
     PlotMarkers -> Table[{mk[[i]], 7}, {i, Length[as]}], PlotStyle -> cs],
    FrameLabel -> {"t", "y(t)"}, PlotRange -> All, Evaluate[frameOpts]]];
f1a = Legended[relaxPanel[{0.3, 0.5, 0.7, 0.9}, cols[[1 ;; 4]], 10., 400],
   mleg[cols[[1 ;; 4]], mk[[1 ;; 4]], ("\[Alpha] = " <> ToString[#]) & /@ {0.3, 0.5, 0.7, 0.9}]];
f1b = Legended[relaxPanel[{1.2, 1.5, 1.8}, cols[[5 ;; 7]], 10., 400],
   mleg[cols[[5 ;; 7]], mk[[1 ;; 3]], ("\[Alpha] = " <> ToString[#]) & /@ {1.2, 1.5, 1.8}]];
save["fig1_mittag_leffler", row[lbl[f1a, "(a) 0 < \[Alpha] < 1, y(0) = 1"], lbl[f1b, "(b) 1 < \[Alpha] < 2, y(0) = 1, y'(0) = 0"]]];

(* ---------------- Fig. 2 and Table: convergence on the Diethelm-Ford-Freed benchmark ---------------- *)
ns = {20, 40, 80, 160, 320, 640};
methods = {"PECE", "PECE3", "TRAP", "LUB1", "LUB2"};
methodNames = <|"PECE" -> "PECE (M = 1)", "PECE3" -> "PECE (M = 3)", "TRAP" -> "implicit trapezoidal",
   "LUB1" -> "Lubich BDF1", "LUB2" -> "Lubich BDF2 + starting weights"|>;
solveWith[m_, f_, a_, y0_, T_, n_] := Switch[m,
   "PECE", FOLPECE[f, a, y0, {T, n}],
   "PECE3", FOLPECE[f, a, y0, {T, n}, "CorrectorIterations" -> 3],
   "TRAP", FOLTrapezoidal[f, a, y0, {T, n}],
   "LUB1", FOLLubich[f, a, y0, {T, n}, "Order" -> 1],
   "LUB2", FOLLubich[f, a, y0, {T, n}]];
benchErr[m_, a_, n_] := Abs[Last[solveWith[m, bench[a], a, y0B[a], 1., n][[2]]] - benchExact[a, 1.]];
alphasTab = {0.25, 0.5, 0.75, 1.25, 1.5, 1.75};
errTab = Association@Flatten@Table[{m, a} -> Table[benchErr[m, a, n], {n, ns}], {m, methods}, {a, alphasTab}];
log["benchmark errors done"];
Export[FileNameJoin[{$tab, "benchmark_errors.csv"}],
  Prepend[Flatten[Table[Join[{m, a}, errTab[{m, a}]], {m, methods}, {a, alphasTab}], 1], Join[{"method", "alpha"}, ToString /@ ns]], "CSV"];
convPanel[a_] := Show[
   ListLogLogPlot[Table[Transpose[{ns, errTab[{m, a}]}], {m, methods}], Joined -> True,
    PlotMarkers -> Table[{mk[[i]], 7}, {i, Length[methods]}], PlotStyle -> cols[[1 ;; Length[methods]]]],
   ListLogLogPlot[{slope[{40, 640}, 3 errTab[{"PECE3", a}][[2]], 2], slope[{40, 640}, 3 errTab[{"LUB1", a}][[2]], 1]},
    Joined -> True, PlotStyle -> Directive[Gray, Dashed]],
   FrameLabel -> {"number of steps N", "|y_N - y(1)|"}, PlotRange -> All, Evaluate[frameOpts]];
f2a = Legended[convPanel[0.5], mleg[cols[[1 ;; 5]], mk[[1 ;; 5]], methodNames /@ methods]];
f2b = Legended[convPanel[1.5], mleg[cols[[1 ;; 5]], mk[[1 ;; 5]], methodNames /@ methods]];
save["fig2_convergence", row[lbl[f2a, "(a) benchmark, \[Alpha] = 0.5 (dashed: slopes 2 and 1)"], lbl[f2b, "(b) benchmark, \[Alpha] = 1.5"]]];
orders = Association@Table[k -> Log2[errTab[k][[-2]]/errTab[k][[-1]]], {k, Keys[errTab]}];
fmt2[x_] := ToString[NumberForm[x, {3, 2}]];
fmtE[x_] := Module[{e = Floor[Log10[x]]}, "$" <> ToString[NumberForm[N[x/10^e], {2, 1}]] <> "\\times10^{" <> ToString[e] <> "}$"];
Export[FileNameJoin[{$tab, "benchmark_orders.tex"}], StringRiffle[Table[
     StringRiffle[Join[{ToString[a]}, Flatten@Table[{fmtE[errTab[{m, a}][[-1]]], fmt2[orders[{m, a}]]}, {m, methods}]], " & "] <> " \\\\",
     {a, alphasTab}], "\n"], "Text"];

(* ---------------- Fig. 3: Caputo-derivative operators ---------------- *)
nsD = {10, 20, 40, 80, 160, 320, 640};
dExp[a_] := N[Sum[1/Gamma[k + 1 - a], {k, Ceiling[a], 80}], 20];
errSF[a_, u_, d0_, ex_] := Table[Abs[Last[FOLCaputoStrongFDM[u, a, 1./n, n, "InitialDerivative" -> d0]] - ex], {n, nsD}];
errL1[a_, p_] := Table[Abs[Last[FOLCaputoL1[N[Range[0, n]/n]^p, a, 1./n]] - FOLCaputoPower[p, a, 1.]], {n, nsD}];
errL1exp[a_] := Table[Abs[Last[FOLCaputoL1[Exp[N[Range[0, n]/n]], a, 1./n]] - dExp[a]], {n, nsD}];
e3a = {errSF[0.5, Exp, 1., dExp[0.5]], errL1exp[0.5], errSF[0.5, #^3.5 &, 0., FOLCaputoPower[3.5, 0.5, 1.]], errL1[0.5, 3.5]};
f3a = Legended[Show[ListLogLogPlot[Table[Transpose[{nsD, e}], {e, e3a}], Joined -> True,
     PlotMarkers -> Table[{mk[[i]], 7}, {i, 4}], PlotStyle -> cols[[{1, 2, 3, 4}]]],
    ListLogLogPlot[{slope[{20, 640}, 2 e3a[[1, 2]], 2], slope[{20, 640}, 2 e3a[[2, 2]], 1.5]}, Joined -> True, PlotStyle -> Directive[Gray, Dashed]],
    FrameLabel -> {"number of steps N", "error at t = 1"}, PlotRange -> All, Evaluate[frameOpts]],
   mleg[cols[[1 ;; 4]], mk[[1 ;; 4]], {"Strong-FDM, e^t", "L1, e^t", "Strong-FDM, t^3.5", "L1, t^3.5"}]];
alphasD = {0.3, 0.8, 1.3, 1.7};
e3b = Table[errSF[a, Exp, 1., dExp[a]], {a, alphasD}];
f3b = Legended[Show[ListLogLogPlot[Table[Transpose[{nsD, e}], {e, e3b}], Joined -> True,
     PlotMarkers -> Table[{mk[[i]], 7}, {i, 4}], PlotStyle -> cols[[{1, 3, 5, 4}]]],
    ListLogLogPlot[{slope[{20, 640}, 2 Max[e3b[[All, 2]]], 2]}, Joined -> True, PlotStyle -> Directive[Gray, Dashed]],
    FrameLabel -> {"number of steps N", "error at t = 1"}, PlotRange -> All, Evaluate[frameOpts]],
   mleg[cols[[{1, 3, 5, 4}]], mk[[1 ;; 4]], ("\[Alpha] = " <> ToString[#]) & /@ alphasD]];
save["fig3_caputo_operators", row[lbl[f3a, "(a) \[Alpha] = 0.5 (dashed: slopes 2 and 2 - \[Alpha])"], lbl[f3b, "(b) Strong-FDM on e^t, 0 < \[Alpha] < 2"]]];
(* orders between N = 160 and 320, before the round-off floor of the h^-p stencil *)
sfOrders = Table[Log2[e[[-3]]/e[[-2]]], {e, e3b}];

(* ---------------- Fig. 4: applications (commensurate systems) ---------------- *)
duff = {#2[[2]], -0.15 #2[[2]] + #2[[1]] - #2[[1]]^3 + 0.3 Cos[#1]} &;
duffSols = Table[FOLPECE[duff, a, {0.5, 0.}, {100., 8000}, "CorrectorIterations" -> 2], {a, {0.9, 0.95, 1.}}];
log["Duffing done"];
f4a = Legended[ListLinePlot[Table[s[[2]][[4001 ;;]], {s, duffSols}], PlotStyle -> (Directive[#, AbsoluteThickness[0.6]] & /@ cols[[{1, 2, 3}]]),
    FrameLabel -> {"x", "y"}, PlotRange -> All, Evaluate[frameOpts]],
   leg[cols[[{1, 2, 3}]], {"\[Alpha] = 0.9", "\[Alpha] = 0.95", "\[Alpha] = 1"}]];
lorenz = {10. (#2[[2]] - #2[[1]]), #2[[1]] (28. - #2[[3]]) - #2[[2]], #2[[1]] #2[[2]] - 8./3 #2[[3]]} &;
lor = FOLPECE[lorenz, 0.995, {1., 1., 1.}, {60., 12000}, "CorrectorIterations" -> 2];
log["Lorenz done"];
f4b = ListLinePlot[lor[[2]][[2001 ;;, {1, 3}]], PlotStyle -> Directive[cols[[1]], AbsoluteThickness[0.4]],
   FrameLabel -> {"x", "z"}, PlotRange -> All, Evaluate[frameOpts]];
save["fig4_applications", row[lbl[f4a, "(a) forced Duffing system, t \[Element] [50, 100]"], lbl[f4b, "(b) Lorenz system, \[Alpha] = 0.995, t \[Element] [10, 60]"]]];

(* ---------------- Fig. 5: variable order (A. H. Salas logistic example) ---------------- *)
alFun = 0.75 + 0.2/(1 + Exp[-1.5 (# - 3)]) &;
logi = 2 #2 (1 - #2) &;
vo = FOLVariableOrderPECE[logi, alFun, 0.2, {10., 2000}];
c1 = FOLPECE[logi, 0.75, 0.2, {10., 2000}]; c2 = FOLPECE[logi, 0.95, 0.2, {10., 2000}];
f5a = Plot[alFun[t], {t, 0, 10}, PlotStyle -> Directive[cols[[4]], AbsoluteThickness[1.6]], FrameLabel -> {"t", "\[Alpha](t)"},
   PlotRange -> {0.7, 1}, Evaluate[frameOpts]];
f5b = Legended[ListLinePlot[{Transpose[vo], Transpose[c1], Transpose[c2]},
    PlotStyle -> {Directive[cols[[4]], AbsoluteThickness[1.8]], Directive[cols[[1]], Dashed], Directive[cols[[2]], Dashed]},
    FrameLabel -> {"t", "y(t)"}, PlotRange -> All, Evaluate[frameOpts]],
   leg[{cols[[4]], Directive[cols[[1]], Dashed], Directive[cols[[2]], Dashed]}, {"variable order \[Alpha](t)", "\[Alpha] = 0.75", "\[Alpha] = 0.95"}]];
save["fig5_variable_order", row[lbl[f5a, "(a) logistic order \[Alpha](t)"], lbl[f5b, "(b) fractional logistic equation"]]];
voRef = Last[FOLVariableOrderPECE[logi, alFun, 0.2, {10., 8000}][[2]]];
voErr = Table[Abs[Last[FOLVariableOrderPECE[logi, alFun, 0.2, {10., n}][[2]]] - voRef], {n, {250, 500, 1000, 2000}}];
voOrder = Log2[voErr[[-2]]/voErr[[-1]]];

(* ---------------- Fig. 6: residual certificate and cost ---------------- *)
aR = 0.5; fR = Function[{t, y}, Gamma[4]/Gamma[4 - aR] t^(3 - aR) + 1 + t^3 - y];
resData = Table[Module[{sol = FOLPECE[fR, aR, 1., {1., n}, "CorrectorIterations" -> 3]},
    {n, Transpose[{Rest[sol[[1]]], Abs[FOLResidual[fR, aR, sol, "InitialDerivative" -> 0.]]}],
     Transpose[{Rest[sol[[1]]], Abs[FOLResidual[fR, aR, sol, "Operator" -> "L1"]]}],
     Max[Abs[sol[[2]] - (1 + sol[[1]]^3)]]}], {n, {40, 160, 640}}];
f6a = Legended[ListLogPlot[Join[resData[[All, 2]], resData[[All, 3]]], Joined -> True,
    PlotStyle -> Join[Directive[#, AbsoluteThickness[1.3]] & /@ cols[[{1, 2, 3}]], Directive[#, Dashed] & /@ cols[[{1, 2, 3}]]],
    FrameLabel -> {"t", "|D^\[Alpha] y_h - f(t, y_h)|"}, PlotRange -> All, Evaluate[frameOpts]],
   leg[cols[[{1, 2, 3}]], ("N = " <> ToString[#]) & /@ {40, 160, 640}]];
nsT = {100, 200, 400, 800, 1600};
timing = Table[{n, First@AbsoluteTiming[solveWith[m, -#2 &, 0.5, 1., 2., n]]}, {m, {"PECE", "TRAP", "LUB2"}}, {n, nsT}];
f6b = Legended[ListLogLogPlot[timing, Joined -> True, PlotMarkers -> Table[{mk[[i]], 7}, {i, 3}], PlotStyle -> cols[[{1, 3, 5}]],
    FrameLabel -> {"number of steps N", "CPU time (s)"}, PlotRange -> All, Evaluate[frameOpts]],
   mleg[cols[[{1, 3, 5}]], mk[[1 ;; 3]], {"PECE", "implicit trapezoidal", "Lubich BDF2"}]];
save["fig6_residual_and_cost", row[lbl[f6a, "(a) residual: Strong-FDM (solid), L1 (dashed)"], lbl[f6b, "(b) run time, scalar problem"]]];
log["figures done"];

(* ---------------- numbers for the manuscript ---------------- *)
texM[x_] := Module[{e = Floor[Log10[N[x]]]}, "\\ensuremath{" <> ToString[NumberForm[N[x/10^e], {2, 1}]] <> "\\times10^{" <> ToString[e] <> "}}"];
testRows = Quiet@Check[Rest@Import[FileNameJoin[{$dir, "test_report.csv"}], "CSV"], {}];
numbers = {
   "wlTests" -> ToString[Length[testRows]], "wlPassed" -> ToString[Count[testRows[[All, 2]], "Success"]],
   "pythonTests" -> "20",
   "mmaVersion" -> StringTrim[First@StringSplit[$Version, "("]],
   "sfOrderMin" -> fmt2[Min[sfOrders]], "sfOrderMax" -> fmt2[Max[sfOrders]],
   "voOrder" -> fmt2[voOrder],
   "resMaxSF" -> texM[Max[resData[[-1, 2, All, 2]]]], "resMaxLone" -> texM[Max[resData[[-1, 3, All, 2]]]],
   "errMaxRes" -> texM[resData[[-1, 4]]],
   "timePECE" -> ToString[NumberForm[timing[[1, -1, 2]], {3, 2}]], "timeTRAP" -> ToString[NumberForm[timing[[2, -1, 2]], {3, 2}]],
   "timeLUB" -> ToString[NumberForm[timing[[3, -1, 2]], {3, 2}]],
   "ordPECEhalf" -> fmt2[orders[{"PECE", 0.5}]], "ordPECEthreeHalf" -> fmt2[orders[{"PECE3", 1.5}]],
   "ordLUBhalf" -> fmt2[orders[{"LUB2", 0.5}]]};
Export[FileNameJoin[{$tab, "numbers.tex"}], StringRiffle[("\\newcommand{\\" <> #[[1]] <> "}{" <> #[[2]] <> "}") & /@ numbers, "\n"] <> "\n", "Text"];
Export[FileNameJoin[{$tab, "environment.txt"}], StringJoin["Generated by MakeFigures.wls with ", $Version, "\nPackage version ", FOLVersion,
   "\nDate ", DateString[], "\n"], "Text"];
log["numbers ", numbers];
log["ALL DONE"];
