ClearAll["Global`*"]
\[Gamma] = 1.;
\[Omega]1 = 1.01;
\[Omega]2 = 0.99;
A0 = 1.;
\[Omega] = 1.;
numberOfPart = 50; (*number of particles in the simulation*)
\
starttime = 0.05; (*controles how the trails are erased*)
runtime = 35;
inc = .025; (*time increments*)
ran := RandomReal[{-10, 10}];
pr = {{-10, 10}, {-10.5, 10.5}};
ftix = {{-3 \[Pi], -2 \[Pi], -\[Pi], 0, \[Pi], 2 \[Pi], 3 \[Pi]}, 
   Automatic};
ar = 1;
is = 600;
InitialConds = Table[Table[ran, 4], numberOfPart];

DiffSol[x0_, y0_, vx0_, vy0_, t1_, t2_] := 
 DiffSol[x0, y0, vx0, vy0, t1, t2] =
  NDSolve[{
    (x^\[Prime]\[Prime])[
      t] == -2 (\[Gamma] Derivative[1][x][t] - 
        A0 \[Omega]1 Cos[\[Omega]1 x[t]] Sin[t \[Omega]1] Derivative[
          1][y][t] - 
        A0 \[Omega]2 Cos[\[Omega]2 x[t]] Sin[t \[Omega]2] Derivative[
          1][y][t]), (y^\[Prime]\[Prime])[
      t] == -2 (A0 \[Omega]1 Cos[t \[Omega]1] Sin[\[Omega]1 x[t]] + 
        A0 \[Omega]2 Cos[
          t \[Omega]2] Sin[\[Omega]2 x[t]] + \[Omega]^2 y[t] + 
        2 A0 \[Omega]1 Cos[\[Omega]1 x[t]] Sin[
          t \[Omega]1] Derivative[1][x][t] + 
        2 A0 \[Omega]2 Cos[\[Omega]2 x[t]] Sin[
          t \[Omega]2] Derivative[1][x][t] + \[Gamma] Derivative[1][
          y][t]),
    x'[0] == vx0, x[0] == x0, y'[0] == vy0, y[0] == y0},
   {x, y}, {t, t1, t2}]

(*create disks which are colored based on their x and y coordinates*)

diskplot[x0_, y0_, vx0_, vy0_, t1_, t2_] :=
 
 Graphics[{Hue[x0 + y0], 
   Disk[Evaluate[{x[t2], y[t2]} /. 
       DiffSol[x0, y0, vx0, vy0, t1, t2]][[1]], .25]}]

makePlot[x0_, y0_, vx0_, vy0_, t1_, t2_] :=
 
 ParametricPlot[
  Evaluate[{x[t], y[t]} /. DiffSol[x0, y0, vx0, vy0, t1, t2]], {t, t1,
    t2},
  PlotRange -> pr,
  AspectRatio -> ar,
  PlotTheme -> "Marketing",
  ImageSize -> is,
  PlotStyle -> Hue[x0 + y0],
  FrameLabel -> {"x", "y"},
  LabelStyle -> {FontFamily -> "Latex", FontSize -> 25},
  FrameTicks -> ftix
  ]

bothPlot2[t1_, t2_] := bothPlot2[t1, t2] =
  Show[
   Table[{
     makePlot[
      InitialConds[[i, 1]],
      InitialConds[[i, 2]],
      InitialConds[[i, 3]],
      InitialConds[[i, 4]],
      t1, t2],
     diskplot[
      InitialConds[[i, 1]],
      InitialConds[[i, 2]],
      InitialConds[[i, 3]],
      InitialConds[[i, 4]],
      t1, t2]},
    {i, numberOfPart}],
   PlotRange -> pr,
   AspectRatio -> ar,
   PlotTheme -> "Scientific",
   ImageSize -> is,
   FrameLabel -> {"x", "y"},
   LabelStyle -> {FontFamily -> "Latex", FontSize -> 25},
   FrameTicks -> ftix]

tableOfPlots = 
  ParallelTable[bothPlot2[0, xx], {xx, starttime, runtime, inc}];

ListAnimate[tableOfPlots]

exportFileName = 
  NotebookDirectory[] <> ToString[numberOfPart] <> "_runtime" <> 
   ToString[runtime] <> "_particles.gif";
Export[exportFileName, tableOfPlots];
