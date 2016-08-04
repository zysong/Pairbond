(* ::Package:: *)

(*
Simplex3Draw` Package
*)
(*
Zhiyuan Song
*)
BeginPackage["Simplex3Draw`"]

Simplex3Draw::usage = "Simplex3Draw[R,S,T,P,q,c,n] draw a stream plot and 
a vector field for the evolutionary dynamics of 3 phenotypes, namely, 
unselective defector, selective cooperator, and selective defector in a 
population. Each player can choose to stay or switch after a round of game. 
The stage game payoff structure is featured by {R, S, T, P}. q is the 
probability that a pair of co-players part from each other in the next 
round by accident. c is the cost to search and form a new pair. n controls 
the resolution of the simplex graph."

Begin["`Private`"]
Simplex3Draw[R_,S_,T_,P_,q_,c_,n_]:=
(*
R,S,T,P represent the payoffs in a stage game: p (CC)=R, p (CD)=S, p (DC)=T, p (DD)=P. q is the probability that a pair of co-players part from each other in the next round by accident. c is the cost to search and form a new pair. n controls the resolution of the simplex graph.
*)
(coords[{x_,y_,z_}]:={x/2+y,x Tan[Pi/3]/2};
(* This is the function to project 3d data to a simplex on 2d triangle space. *)
xlist={};
(* Set up an empty array *)
For[i=1,i<=n,i++,
For[j=1;x=i/n,j<n-i,j++,y=j/n;z=1-x-y;xlist=Append[xlist,{x,y,z}]
]
];
(* Uniformly distributed points in a 3D space limited by x+y+z=1, x>=0, y>=0, z>=0. *)
frame=Graphics[{Thick,Line[{{0,0},{1,0},{0.5,Sqrt[3]/2},{0,0}}]},ImageSize->300];
(* The frame of a simplex *)
vertices=Graphics[{Text[Style["Selective Defector",Blue,Medium,Bold],{0-0.02,0-0.02}],Text[Style["Selective Cooperator",Green,Medium,Bold],{1+0.02,0-0.02}],Text[Style["Nonselective Defector",Red,Medium,Bold],{0.5,Sqrt[3]/2+0.02}]},ImageSize->300];
(* Label the three vertices of a simplex *)

seqns={s2 (s1+s2+q s3+q s4)==q (s1+s2+s3+s4)x2,
s3 (s1+q s2+s3+q s4)==q (s1+s2+s3+s4)x3,
s4 (s1+q (s2+s3+s4))==q (s1+s2+s3+s4)x4,
s1==q x1};
(* {s1, s2, s3, s4} is the stationary fractions of "single" players from each type of player pool, namely, unselective cooperator, unselective defector, selective cooperator, selective defector, given the fraction of each type of player in the whole population as {x1, x2, x3, x4}.
*)
u1=(R (s1+s3)+S (s2+s4))/(s1+s2+s3+s4)-q c;
u2=(T (s1+q s3)+P (s2+q s4)-q c (s1+s2+s3+s4))/(s1+s2+q s3+q s4);
u3=(R (s1+ s3)+q S (s2+s4)-q c (s1+s2+s3+s4))/(s1+q s2+s3+q s4);
u4=(T (s1+q s3)+q P (s2+s4)-q c (s1+s2+s3+s4))/(s1+q s2+q s3+q s4);
(* The mean fitness of each type of player at equilibrium *)

eqmx2=Solve[{seqns[[1]],seqns[[2]],u2==u3}/.{x3->1-x2,s1->0,s4->0},{x2,s2,s3}];
eqmx4=Solve[{seqns[[3]],seqns[[2]],u4==u3}/.{x3->1-x4,s1->0,s2->0},{x4,s4,s3}];
stableeqm=coords[#]& /@{{1,0,0},{0,1-x4/.eqmx4[[1]][[1]],x4/.eqmx4[[1]][[1]]}};
unstableeqm=coords[#]& /@{{0,1,0},{0,0,1},{x2/.eqmx2[[1]][[1]],1-x2/.eqmx2[[1]][[1]],0},{x2/.eqmx2[[2]][[1]],1-x2/.eqmx2[[2]][[1]],0},{0,1-x4/.eqmx4[[2]][[1]],x4/.eqmx4[[2]][[1]]}};
points=Graphics[{Circle[#,.01]& /@unstableeqm,Disk[#,.01]& /@stableeqm},ImageSize->300];
(* Label the stable (disk) and unstable (circle) steady states on the edges and vertices. But note that these only work for cases when c > 0. *)

ulist={};
For[i=1,i<=Length[xlist],i++,
uvector={u2,u3,u4}/.Last[Solve[seqns/.{x1->0,x2->xlist[[i]][[1]],x3->xlist[[i]][[2]],x4->xlist[[i]][[3]]}, {s1,s2,s3,s4}]];
ulist=Append[ulist,uvector]];
(* Compute the fitness of each type of player except type 1 as we assume x1 = 0 here. *)

ubarlist=Total[Transpose[ulist*xlist]];
(* Mean fitness of the three types of players (2, 3, 4) *)
gradient3d=xlist*(ulist-ubarlist);
(* Compute the gradient at each point in a 3D space according to the replicator equation. *)

simdata=coords[#]& /@xlist;
(* Project the points from 3D to a simplex *)
simvector=coords[#]& /@gradient3d;
(* Project the gradients from 3D to a simplex *)
vector={};
For[
i=1,i<=Length[simdata],i++,vector=Append[vector,{simdata[[i]],simvector[[i]]}]
];
(* Group each point with the corresponding gradient *)

vectorfield=ListVectorPlot[vector,VectorPoints->simdata,Frame->False];
streamplot=ListStreamPlot[vector,StreamPoints->simdata,Frame->False,RegionFunction->Function[{x,y},Sqrt[3]*x+y<Sqrt[3] && y<Sqrt[3]*x]];
showstreamplot=Show[frame,vertices,streamplot,points, ImageSize->600]
)

End[]
EndPackage[]
