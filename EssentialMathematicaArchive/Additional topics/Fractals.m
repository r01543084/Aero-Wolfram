BeginPackage["Fractals`"];

(* Load required packages *)

Needs["Utilities`FilterOptions`"];

(* Usage messages for exported symbols *)

mandelbrotRate::usage = "mandelbrotRate[cx,cy,lim] computes the escape rate \
for c = cx + I cy under the quadratic map z\[Rule]\!\(z\^2\)+c beginning with \
z=0.  The maximum number of iterations is given by lim, which must be an \
integer."; 

juliaRate::usage = "juliaRate[zx,zy,cx,cy,lim] computes the escape rate for z \
= zx + I zy under the quadratic map z\[Rule]\!\(z\^2\)+c where c = cx +I cy \
is the complex parameter characterizing the map.  The maximum number of \
iterations is given by lim, which must be an integer."; 

DisplayKoch::usage = "DisplayKoch[polygon,n,opts] produces a Koch curve of \
order n based upon the list of points in polygon. By default, the triangle \
rule is applied in the outward direction assuming that the polygon is defined \
counterclockwise. \
\n  Options: \
\n    KochType\[Rule]\"sieve\" applies the triangle-building rule in reverse. \ 
\n    ShowIterations\[Rule]True can be used to display intermediate iterations. \
\n    CloseFigure\[Rule]True connects the last point to the first to close the \
figure."; 

ShowIterations::usage = "a True|False option for DisplayKoch."

KochType::usage = "an option for DisplayKoch: \
\n   \"snowflake\" produces outward triangles \
\n   \"sieve\" produces inward triangles"

CloseFigure::usage = "a True|False option for DisplayKoch that can be used \
when the list of points is not a closed polygon."


Begin["`Private`"];

(* Definitions *)

mandelbrotRate = 
   Compile[{cx, cy, {lim, _Integer}}, 
    Module[{c, z, ct = 0}, c = cx + I*cy; z = c; 
      While[Abs[z] < 2. && ct <= lim, ++ct; z = z*z + c; ]; ct]]; 

juliaRate = Compile[{zx, zy, cx, cy, {lim, _Integer}}, 
    Module[{c, z, ct = 0}, c = cx + I*cy; z = zx + I*zy; 
      While[Abs[z] < 2. && ct <= lim, ++ct; z = z*z + c; ]; ct]]; 

perp[vector:{x_, y_}] := {y, -x}/Sqrt[x^2 + y^2]; 
perp[pt1:{x1_, y1_}, pt2:{x2_, y2_}] := perp[pt1 - pt2]; 

newpt1[pt1:{x1_, y1_}, pt2:{x2_, y2_}] := pt1 + 1/3*(pt2 - pt1); 

newpt2[pt1:{x1_, y1_}, pt2:{x2_, y2_}] := 
   1/2*(pt1 + pt2) + 1/Sqrt[12]*perp[pt1, pt2]*Sqrt[(pt2 - pt1) . (pt2 - pt1)]
 
newpt3[pt1:{x1_, y1_}, pt2:{x2_, y2_}] := pt1 + 2/3*(pt2 - pt1); 

KochRule = Line[{pt1:{x1_, y1_}, pt2:{x2_, y2_}}] :> 
    {Line[{pt1, newpt1[pt1, pt2]}], 
     Line[{newpt1[pt1, pt2], newpt2[pt1, pt2]}], 
     Line[{newpt2[pt1, pt2], newpt3[pt1, pt2]}], 
     Line[{newpt3[pt1, pt2], pt2}]}; 
KochSplit[x_] := x /. KochRule

DisplayKoch[polygon:{{_, _}..}, n_, opts___Rule] := 
Module[{curve, type, showIterations, close, NestFunction, optKoch, optPlot}, 
	If[n<=0,Message[DisplayKoch::negiter];Abort[]];
	{type, showIterations, close} = 
		{KochType, ShowIterations, CloseFigure} /. {opts} /.Options[DisplayKoch]; 
	curve = If[type == "sieve", Reverse[polygon], polygon]; 
	If[close, AppendTo[curve, First[curve]]]; 
	curve = Line /@ Partition[curve, 2, 1]; 
	NestFunction = If[showIterations, NestList, Nest]; 
	optPlot = FilterOptions[Graphics, opts]; 
	Show[Graphics[NestFunction[KochSplit, curve, n]], optPlot, 
		AspectRatio -> Automatic]
]; 

(* Options for exported functions *)

Options[DisplayKoch] = 
   {ShowIterations -> False, KochType -> "snowflake", CloseFigure -> False}; 

(* Error conditions *)

DisplayKoch::negiter="negative iteration number";

DisplayKoch::badarg="`1` arguments given where 2 were expected";

DisplayKoch::notpoly="first argument should be a list of pairs";

DisplayKoch[x_,n_,opts___]:=Message[DisplayKoch::notpoly];

DisplayKoch[x_]:=Message[DisplayKoch::badarg,1];

End[];

(* Protect exported functions *)

Protect[juliaRate,mandelbrotRate,DisplayKoch];

EndPackage[]