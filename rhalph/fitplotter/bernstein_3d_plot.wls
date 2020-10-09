#!/usr/bin/env wolframscript

BernsteinPolys[nRho_, nPt_, pars_, rho_, pt_] := Sum[Sum[pars[[l + 1]][[k + 1]]*BernsteinBasis[nRho, k, rho]*BernsteinBasis[nPt, l, pt], {l, 0, nPt}], {k, 0, nRho}]

TF[nRho_, rho_, nPt_, pt_, pars_] := BernsteinPolys[nRho, nPt, pars, rho, pt]

 PlotTF[nRho_, nPt_, pars_] := 
  Plot3D[
  {TF[nRho, \[Rho], nPt, Subscript[p, T], pars], 1.0},
  {\[Rho], 1, 0}, {Subscript[p, T], 0, 1},
  PlotRange -> {{0, 1}, {0, 1}, {-2, 3}},
  AxesLabel -> {"\!\(\*SubscriptBox[\(m\), \(SD\)]\) [GeV]", "\!\(\*SubscriptBox[\(p\), \(T\)]\) [GeV]","\!\(\*SubscriptBox[\(R\), \(p/f\)]\)"},
  Ticks -> {
        Table[{i, 50 + (200 - 50) i}, {i, 0, 1, 1/7.5}],
        {{0.07142857, 550}, {0.14285714, 600}, {0.25,675}, {0.42857143, 800}, {1., 1200}},
        {-1, 0,1, 2}},
  PlotStyle -> {{Red, Opacity[0.7]}, LightGray},
  LabelStyle -> Small
  ]

argv = Rest @ $ScriptCommandLine;
plot = PlotTF[ToExpression[argv[[2]]],ToExpression[argv[[3]]],ToExpression[argv[[4]]]];

Export[argv[[1]],plot,ImageSize -> 600, "CompressionLevel" -> 0]
