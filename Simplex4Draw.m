(* ::Package:: *)

(*
Simplex4Draw` Package
*)
(*
Zhiyuan Song
Department of Biology
Stanford University
*)
(*
Implementation
*)
(*
Begin Package
*)
BeginPackage["Simplex4Draw`"]

Simplex4Draw::usage = "The Simplex4Draw` package includes functions to project \
the evolutionary dynamics of 4 phenotypes (cooperator, defector, switching \
cooperator, and switching defector) to a 3-D simplex."

(*
Usage Messages
*)
coords43::usage = "coords43[{x1_,x2_,x3_,x4_}] projects a 4D data-point to a \
3D-simplex."

TwoTierEvol::usage = "TwoTierEvol[X0list_,mor_,q_,c_,R_,S_,T_,P_,mu_,tlimit_] \
computes the evolutionary trajectroy from the given initial conditions X0list,\
and the payoff matrix of a stage game is featured by {R, S, T, P}."

S3ListDraw::usage = "S3ListDraw[DataList_] plots the data in a 3D-simplex."

Begin["`Private`"]
coords43[{x1_,x2_,x3_,x4_}]:={(x2+x1)/2+x3,(x2+x1/3) Tan[Pi/3]/2,Sqrt[2/3] x1}

TwoTierEvol[X0list_,mor_,q_,c_,R_,S_,T_,P_,mu_,tlimit_]:=
Module[{xrecord3d={},xlist,
n=Length[X0list],
switch={{q,q,q,q},{q,q,1,1},{q,1,q,1},{q,1,1,1}},
rstp={{R,S},{T,P}},
mutation={{(1-mu)^2, mu (1-mu), mu (1-mu), mu^2}, {mu (1-mu), (1-mu)^2, mu^2, mu (1-mu)}, 
{mu (1-mu), mu^2, (1-mu)^2, mu (1-mu)}, {mu^2, mu (1-mu), mu (1-mu), (1-mu)^2}},
payoff,t,gradientx,wlist,wbar,SingleFrac,single,
xrecord,xrecordlist={},PairMatrix,TransMatrix},
payoff=ArrayFlatten[{{rstp,rstp},{rstp,rstp}}]-c*switch;
For[i=1,i<=n,i++,
xlist=X0list[[i]];PairMatrix=ConstantArray[xlist,4];
xrecord={xlist};
t=0;
While[
t<tlimit,
wlist=Table[payoff[[i]].PairMatrix[[i]],{i,4}];
wbar=xlist.wlist;
gradientx=mor*(Table[Sum[mutation[[i]][[j]]*xlist[[j]]*wlist[[j]]/wbar,{j,4}],{i,4}]-xlist);
SingleFrac=Table[switch[[i]].PairMatrix[[i]],{i,4}];
single=SingleFrac*xlist+gradientx;
xlist=gradientx+xlist;
xrecord=Append[xrecord,xlist];
TransMatrix=Table[Table[switch[[k]][[i]]*single[[j]]/Total[single],{i,4},{j,4}]+
IdentityMatrix[4]-DiagonalMatrix[switch[[k]]],{k,4}];
PairMatrix=Table[PairMatrix[[i]].TransMatrix[[i]],{i,4}];
t++
];
xrecordlist=Append[xrecordlist,xrecord];
];
xrecordlist
]

S3ListDraw[DataList_]:=
Module[{frame,vertices,trajectory,xrecord3d},
frame=Graphics3D[{Thick,Line[{{0,0,0},{1,0,0},{0.5,Sqrt[3]/2, 0},
{0,0, 0},{1/2,Sqrt[3]/6,Sqrt[6]/3},{1,0,0},
{0.5,Sqrt[3]/2, 0},{1/2,Sqrt[3]/6,Sqrt[6]/3}}]}, 
Boxed->False, ImageSize->300];
(* The frame of a simplex *)
vertices=Graphics3D[{Text[Style["SD",Medium,Bold],{0-0.02,0-0.02,0-0.02}],
Text[Style["SC",Medium,Bold],{1+0.02,0-0.02, 0-0.02}], 
Text[Style["ND",Medium,Bold],
{0.5,Sqrt[3]/2+0.02, 0-0.02}], 
Text[Style["NC",Medium,Bold], 
{0.5,Sqrt[3]/6,Sqrt[6]/3+0.02 }]},Boxed->False,ImageSize->300];
(* Label the three vertices of a simplex *)
xrecord3d=Table[coords43[#]& /@DataList[[i]],{i,Length[DataList]}];
trajectory=ListPointPlot3D[xrecord3d,Axes->False,Boxed->False,ImageSize->300];
Show[frame,vertices,trajectory, ImageSize->600]
]

End[]

Protect[Evaluate[Context[] <> "*"]]
EndPackage[]
