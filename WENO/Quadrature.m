(*:Mathematica Version: 6.0 *)

(*:Package Version: 1.02 *)

(*:Name: MathWorld`Quadrature` *)

(*:Author: Eric W. Weisstein *)

(*:URL:
  http://mathworld.wolfram.com/packages/
*)

(*:Summary:
*)

(*:History:
  v1.00 (2002-01-31): Written
  v1.01 (2003-10-19): context updated
  v1.02 (2005-09-15): *Exact removed and GaussianQuadrature rewritten
  
  (c) 2002-2007 Eric W. Weisstein
*)

(*:Keywords:
  
*)

(*:Requirements: None. *)

(*:Discussion:
*)

(*:References: *)

(*:Limitations: None known. *)

BeginPackage["MathWorld`Quadrature`"]

ChebyshevQuadrature::usage =
"ChebyshevQuadrature."

GaussianQuadrature::usage =
"GaussianQuadrature[f,W,{x0,x1},n] gives {x,w}, where x and w are lists of the n abscissas and \
weights for Gaussian quadrature for an orthogonal polynomial f of degree n on the interval {x0,x1} \
with weight W[#]&."

LobattoQuadrature::usage =
"LobattoQuadrature."

RadauQuadrature::usage =
"RadauQuadrature."


Begin["`Private`"]

A[m_,ChebyshevT]:=2^(m-1)
gamma[m_,ChebyshevT]:=Pi/2

A[m_,HermiteH]:=2^m
gamma[m_,HermiteH]:=Sqrt[Pi] 2^m m!

A[m_,LegendreP]:=(2m)!/2^m/(m!)^2
gamma[m_,LegendreP]:=2/(2m+1)

A[m_,LaguerreL]:=(-1)^m/m!
gamma[m_,LaguerreL]:=1

(*A[m_,f]:=CoefficientList[f[m,t],t][[-1]];
  gamma[m_,f]:=NIntegrate[wt[t] (f[m,t])^2,{t,range[[1]],range[[2]]}];*)

G[n_,x_]:=Module[{m,g,y},
    Expand[x^n Normal[
          Series[Exp[(-2+Log[1-y](1-1/y)+Log[1+y](1+1/y))n/2],{y,0,n}]/.{y->1/x}]]]

ChebyshevQuadrature[n_]:=Module[{x},
    x=Roots[G[n,t]==0,t];
    Table[x[[i,2]],{i,1,n}]
]

GaussianQuadrature[f_,wt_,range_List,n_]:=Module[{x,w,t},
    x=Root[f[n,t],t,#]&/@Range[n];
 	w=gamma[n-1,f]A[n,f]/A[n-1,f]/f[n-1,#]/Derivative[0,1][f][n,#]&/@x;
	{x,w}
]

LobattoQuadrature[n_]:=
  Module[{i,w,xor,x},If[n>3,xor=Roots[D[LegendreP[n-1,t],t]==0,t];
      x=Table[xor[[i,2]],{i,1,n-2}],x={0};];
    w=Table[2/n/(n-1)/LegendreP[n-1,x[[i]]]^2,{i,1,n-2}];
    x=Join[{-1},x,{1}];
    w=Join[{2/n/(n-1)},w,{2/n/(n-1)}];
    {x,w//N}]

RadauQuadrature[n_]:=
  Module[{i,w,xor,x},
    xor=Roots[Together[(LegendreP[n-1,t]+LegendreP[n,t])/(1+t)]==0,t];
    x=Table[If[n==2,xor[[2]],xor[[i,2]]],{i,1,n-1}];
    w=Table[(1-x[[i]])/n^2/LegendreP[n-1,x[[i]]]^2,{i,1,n-1}];
    x=Join[{-1},x];
    w=Join[{2/n^2},w];
    {x,w//N}]

End[]

EndPackage[]

(* Protect[  ] *)
