BeginPackage["Curfit`"]; 

(* Load required packages *)

(* Usage messages for exported functions *)

Linfit::usage = "Linfit[data,f,var] performs linear least-squares fitting to \
data with one independent variable \!\(\* \
StyleBox[\"x\",\nFontSlant->\"Italic\"]\) and one dependent variable \!\(\* \
StyleBox[\"y\",\nFontSlant->\"Italic\"]\).  The data array should be a list \
of the form \!\(\* StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"}\",\nFontSlant->\"Italic\"]\), \!\(\* \
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"dy\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"}\",\nFontSlant->\"Italic\"]\), or \!\(\* \
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"dx\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"dy\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"}\",\nFontSlant->\"Italic\"]\) for each measurement.  The list of \
basis functions is \!\(\* StyleBox[\"f\",\nFontSlant->\"Italic\"]\) and \
\!\(\* StyleBox[\"var\",\nFontSlant->\"Italic\"]\) is the symbol used in \
\!\(\* StyleBox[\"f\",\nFontSlant->\"Italic\"]\) for the independent \
variable."; 

Curfit::usage = "Curfit[data,f,var,par] performs nonlinear least-squares \
fitting using the Levenberg-Marquardt algorithm. The data array should be a \
list of the form \!\(\* StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"}\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"dy\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"}\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"or\",\nFontSlant->\"Italic\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"{\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"dx\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"dy\",\nFontSlant->\"Italic\"]\)\!\(\* \
StyleBox[\"}\",\nFontSlant->\"Italic\"]\) for each measurement. The fitting \
function \!\(\* StyleBox[\"f\",\nFontSlant->\"Italic\"]\) is expressed in \
terms of the independent variable \!\(\* \
StyleBox[\"var\",\nFontSlant->\"Italic\"]\) and the parameters \!\(\* \
StyleBox[\"par\",\nFontSlant->\"Italic\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"are\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"given\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"as\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"a\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"list\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"containing\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"{\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"name\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\",\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"initial\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"value\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"}\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"for\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"each\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"adjustable\",\nFontSlant->\"Plain\"]\)\!\(\* StyleBox[\" \
\",\nFontSlant->\"Plain\"]\)\!\(\* \
StyleBox[\"parameter\",\nFontSlant->\"Plain\"]\)."; 

NitsDx::usage = 
   "Linfit or Curfit option for maximum number of iterations wrt dx."; 

ConvDx::usage = "Linfit or Curfit option specifying relative change in \
chisquare for convergence wrt dx."; 

MaxNITS::usage = "Curfit option for maximum number of iterations."; 

Conv::usage = "Curfit option gives maximum change in reduced chisquare \
for convergence."; 

PrintIterations::usage = "Curfit option for printing iterations."; 

IncludeSecondDerivatives::usage = "Curfit option governing use of second \
derivatives in curvature matrix.";

NITS::usage = "number of iterations";
DxNits::usage = "iteration history with uncertainties in independent variable";
Parameters::usage = "best-fit parameters";
Uncertainties::usage = "parameter uncertainties";
StandardDeviation::usage = "standard deviation";
ChiSquare::usage = "final chisquare";
ReducedChiSquare::usage = "reduced chisquare";
Correlations::usage = "correlation matrix for Curfit" ;

Begin["`Private`"]; 

(* Definitions *)

Linfit[data:{{_, _}..} | {{_, _, _}..}, (f_)?VectorQ, var_Symbol, 
    opts___Rule] := linfit[data, f, var]; 

Linfit[data:{{_, _, _, _}..}, (f_)?VectorQ, var_Symbol, opts___Rule] := 
  Module[{x, dx, y, dy, \[Sigma], ndata, ndim, nparam, i, j, params, slopes, fp, 
    errors, yfit, redchi, oldparams, nits, mnits, change, conv}, 
   {ndata, ndim} = Dimensions[data]; If[ndata <= Length[f] || ndim != 4, 
     Message[Linfit::BadDim, {ndata, ndim}]; Return[]]; 
    {x, dx, y, dy} = Transpose[data]; nparam = Length[f]; 
    slopes = Table[0, {ndata}]; fp = D[f, var]; 
    {mnits, conv} = {NitsDx, ConvDx} /. {opts} /. Options[Linfit]; 
    oldparams = Table[0, {i, nparam}]; For[nits = 0; change = 1, 
     nits < mnits && change > conv, nits++, 
     \[Sigma] = Sqrt[dy^2 + (slopes*dx)^2]; {params, errors} = 
       {Parameters, Uncertainties} /. linfit[Transpose[{x, y, \[Sigma]}], f, 
         var]; yfit = params . f /. var -> x; 
      slopes = params . fp /. var -> x; 
      change = Max[Abs[(params - oldparams)/params]]; oldparams = params]; 
    If[mnits > 1 && change > conv, Message[Linfit::conv, 
      NumberForm[change, 2]]]; redchi = \[Sigma]^(-2) . (y - yfit)^2/
      (ndata - nparam); {Parameters -> params, Uncertainties -> errors, 
     ReducedChiSquare -> redchi, NITS -> nits}]

linfit[data_, (f_)?VectorQ, var_] := 
  Module[{x, y, dy, w, M, b, \[CurlyEpsilon], basis, ndata, nparam, ndim, i, j, params, 
    yfit, redchi, errors, stdev}, {ndata, ndim} = Dimensions[data]; 
    If[ndata <= Length[f] || ndim < 2 || ndim > 3, 
     Message[Linfit::BadDim, {ndata, ndim}]; Return[]]; 
    If[ndim == 2, {x, y} = Transpose[data]; w = Table[1, {ndata}], 
     {x, y, dy} = Transpose[data]; w = 1/dy^2]; nparam = Length[f]; 
    M = Table[0, {i, nparam}, {j, nparam}]; b = Table[0, {i, nparam}]; 
    Do[basis = f /. var -> x[[i]]; 
      M = M + w[[i]]*Outer[Times, basis, basis]; 
      b = b + w[[i]]*y[[i]]*basis, {i, 1, ndata}]; \[CurlyEpsilon] = PseudoInverse[M]; 
    params = \[CurlyEpsilon] . b; errors = Table[Sqrt[\[CurlyEpsilon][[i,i]]], {i, nparam}]; 
    yfit = params . f /. var -> x; redchi = w . (y - yfit)^2/
      (ndata - nparam); stdev = Sqrt[redchi]; 
    If[ndim == 2, {Parameters -> params, Uncertainties -> errors*stdev, 
      StandardDeviation -> stdev}, {Parameters -> params, 
      Uncertainties -> errors, ReducedChiSquare -> redchi}]]

Curfit[data:{{_, _}..} | {{_, _, _}..}, f_, var_Symbol, 
   par:{{_, _?NumericQ}..}, opts___Rule] := curfit[data, f, var, par, opts]

Curfit[data:{{_, _, _, _}..}, f_, var_Symbol, pars:{{_, _?NumericQ}..}, 
   opts___Rule] := Module[{x, dx, y, dy, \[Sigma], ndata, ndim, nparam, i, j, 
    par, start, params, yfit, slopes, fp, errors, redchi, oldparams, 
    newparams, nits, mnits, printIterations, dxnits, change, conv, 
    result}, {ndata, ndim} = Dimensions[data]; 
    If[ndata <= Length[par] || ndim != 4, 
     Message[Curfit::BadDim, {ndata, ndim}]; Return[]]; 
    {x, dx, y, dy} = Transpose[data]; {par, start} = Transpose[pars]; 
    oldparams = start; nparam = Length[pars]; 
    slopes = Table[0, {ndata}]; fp = D[f, var]; 
    {mnits, conv, printIterations} = 
     {NitsDx, ConvDx, PrintIterations} /. {opts} /. Options[Curfit]; 
    dxnits = {}; For[nits = 0; change = 1, nits < mnits && 
      change > conv, nits++, \[Sigma] = Sqrt[dy^2 + (slopes*dx)^2]; 
      result = curfit[Transpose[{x, y, \[Sigma]}], f, var, 
        Transpose[{par, oldparams}], opts]; {params, errors} = 
       {Parameters, Uncertainties} /. result; AppendTo[dxnits, 
       NITS /. result]; If[printIterations, Print[params]]; 
      yfit = f /. params /. var -> x; slopes = fp /. params /. var -> x; 
      newparams = par /. params; change = 
       Max[Abs[(newparams - oldparams)/newparams]]; 
      oldparams = newparams]; If[mnits > 1 && change > conv, 
     Message[Curfit::dx, NumberForm[change, 2]]]; 
    redchi = \[Sigma]^(-2) . (y - yfit)^2/(ndata - nparam); 
    AppendTo[result, DxNits -> dxnits]]

curfit[data_, f_, var_, pars_, opts___Rule] := 
  Module[{par, start, npar, ndata, ndim, nits, a, da, \[Lambda], x, y, dy, w, 
    yfit, derivs, dmat, deriv2, dmat2, \[Beta], \[Alpha], \[CurlyEpsilon], errors, corr, chisq, 
    redchi, stdev, change, conv, mnits, i, j, oldchi, newchi, 
    printIterations, includeSecondDerivatives}, 
   npar = Length[pars]; {par, start} = Transpose[pars]; 
    {ndata, ndim} = Dimensions[data]; If[ndata <= npar || ndim < 2 || 
      ndim > 3, Message[Curfit::BadDim, {ndata, ndim}]; Return[]]; 
    If[ndim == 2, {x, y} = Transpose[data]; w = Table[1, {ndata}], 
     {x, y, dy} = Transpose[data]; w = 1/dy^2]; 
    {mnits, conv, printIterations, includeSecondDerivatives} = 
     {MaxNITS, Conv, PrintIterations, IncludeSecondDerivatives} /. {opts} /. 
      Options[Curfit]; derivs = Table[Simplify[D[f, par[[i]]]], {i, npar}]; 
    If[includeSecondDerivatives, deriv2 = 
      Table[Simplify[D[f, par[[i]], par[[j]]]], {i, npar}, {j, npar}]]; 
    change = 1.; \[Lambda] = 0.001; For[nits = 0; a = start, 
     nits <= mnits && change > conv, nits++, 
     yfit = f /. Thread[par -> a] /. var -> x; oldchi = w . (y - yfit)^2; 
      dmat = Table[derivs /. Thread[par -> a] /. var -> x[[i]], {i, ndata}]; 
      \[Beta] = (w*(y - yfit)) . dmat; \[Alpha] = Transpose[w*dmat] . dmat; 
      If[includeSecondDerivatives, 
       dmat2 = Table[deriv2 /. Thread[par -> a] /. var -> x[[i]], 
          {i, ndata}]; \[Alpha] = \[Alpha] - (w*(y - yfit)) . dmat2]; 
      Do[\[Alpha][[i,i]] = (1 + \[Lambda])*\[Alpha][[i,i]], {i, npar}]; 
      \[CurlyEpsilon] = PseudoInverse[\[Alpha]]; da = \[CurlyEpsilon] . \[Beta]; 
      yfit = f /. Thread[par -> a + da] /. var -> x; 
      newchi = w . (y - yfit)^2; If[newchi < oldchi, 
       \[Lambda] = \[Lambda]/10; a = a + da; change = (oldchi - newchi)/oldchi, 
       \[Lambda] = 10*\[Lambda]]; If[printIterations, Print[{oldchi, newchi, \[Lambda], 
         change, a + da}]]; ]; If[mnits > 1 && change > conv, 
     Message[Curfit::conv, NumberForm[change, 2], NumberForm[\[Lambda], 2]]]; 
    \[Alpha] = Transpose[w*dmat] . dmat; If[includeSecondDerivatives, 
     dmat2 = Table[deriv2 /. Thread[par -> a] /. var -> x[[i]], {i, ndata}]; 
      \[Alpha] = \[Alpha]*(w*(y - yfit)) . dmat2]; \[CurlyEpsilon] = PseudoInverse[\[Alpha]]; 
    errors = Table[Sqrt[\[CurlyEpsilon][[i,i]]], {i, npar}]; 
    corr = Table[\[CurlyEpsilon][[i,j]]/(errors[[i]]*errors[[j]]), {i, npar}, 
      {j, npar}]; yfit = f /. Thread[par -> a] /. var -> x; 
    chisq = w . (y - yfit)^2; redchi = chisq/(ndata - npar); 
    stdev = Sqrt[redchi]; If[ndim == 2, {Parameters -> Thread[par -> a], 
      Uncertainties -> errors*stdev, StandardDeviation -> stdev}, 
     {Parameters -> Thread[par -> a], Uncertainties -> errors, 
      ReducedChiSquare -> redchi, ChiSquare -> chisq, NITS -> nits, 
      Correlations -> corr}]]

(* Options for exported functions *)

Options[Linfit] = {NitsDx -> 10, ConvDx -> 0.001}; 

Options[Curfit] = {NitsDx -> 10, ConvDx -> 0.001, MaxNITS -> 100, 
    Conv -> 0.0001, PrintIterations -> False, IncludeSecondDerivatives -> 
     False}; 

(* Error conditions *)

Linfit::BadDim = "Bad dimensions in data array: `1`"; 

Linfit::conv = 
   "Warning: failed to converge.  Maximum relative change = `1`."; 

Curfit::BadDim = "Bad dimensions in data array: `1`"; 

Curfit::conv = 
   "Warning: failed to converge.  Relative change = `1`; \[Lambda] = `2`"; 

Curfit::dx = "Warning: dx iteration failed to converge.  Maximum \
relative change = `1`"; 

End[]; 

(* Protect exported functions *)

Protect[Linfit, Curfit]; 

EndPackage[]