(* ::Package:: *)

(*
ABS` Package
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
BeginPackage["ABC`"]

ABC::usage= "The ABC` Package includes functions to run agent-based 
simulations of the coevolution of cooperation and partner switching."

ABS::usage= "ABS[n0_,xt0_,mortality_,q_,c_,R_,S_,T_,P_,mu_,step_,run_] 
returns the evolving trajectories of {x11, x21, x12, x22}. 
n0 is the population size; 
xt is the vector of x11, x21, x12 and x22, the frequencies of 'CN', 'DN' 'CS' and 'DS'; 
w is the intensity of selection; 
q is the accidental break-up rate of a pair;
c is the cost of partner switching; 
RSTP are the payoffs of a stage game; 
mu is the mutation rate; 
step is the full length of each run."

Begin["`Private`"]
ABS[n0_,xt0_,mortality_,q_,c_,R_,S_,T_,P_,mu_,step_,run_]:= 
Module[{n=2*Floor[n0/2],n11,n21,n12,n22,partner,traits,trait1,trait2,type,fitness,
single,newpair,payoffmatrix={{R,S},{T,P}},
switch={{0,0,0,0},{0,0,1,1},{0,1,0,1},{0,1,1,1}},
death,partnerdeath,error,partnererror,Issingle,Ndeath,Pdeath,newborn,newtraits,
mutation,xt,xtlist,xlist={}},
If[Total[xt0]>1,Print["Error: x11+x12+x21+x22>1"],
(*Error information*)
Do[
xt=xt0;
xtlist={xt};
n11=Floor[xt[[1]]*n];n21=Floor[xt[[2]]*n];n12=Floor[xt[[3]]*n];n22=n-n11-n21-n12;
traits=Join[ConstantArray[{1,1},n11],ConstantArray[{2,1},n21], ConstantArray[{1,2},n12],ConstantArray[{2,2},n22]];
trait1=Transpose[traits][[1]];
trait2=Transpose[traits][[2]];
type=trait1+trait2*2-2;
(* Code traits 11, 21, 12, 22 as type 1, 2, 3, 4.*)
partner=ConstantArray[0,n];
single=Range[n];
(* Initial conditions for each run *)
Do[
fitness=ConstantArray[0,n];
(* Initial fitness for each step *)
newpair=Partition[RandomSample[single],2];
(* Random matching *)
fitness[[single]]=fitness[[single]]-c;
(* Cost of new pair matching *)
For[
i=1,i<=Length[newpair],i++,
partner[[newpair[[i]][[1]]]]=newpair[[i]][[2]];
partner[[newpair[[i]][[2]]]]=newpair[[i]][[1]]
];
fitness=fitness +Table[payoffmatrix[[trait1[[i]]]][[trait1[[partner[[i]]]]]],{i,n}];
death=RandomInteger[BernoulliDistribution[mortality],n];
(* Random deaths *)
partnerdeath=Table[death[[partner[[i]]]],{i,n}];
error=RandomInteger[BernoulliDistribution[1-Sqrt[1-q]/(1-mortality)],n];
(* Random error *)
partnererror=Table[error[[partner[[i]]]],{i,n}];
Issingle=1-(1-Table[switch[[type[[i]]]][[type[[partner[[i]]]]]],{i,n}])(1-death)(1-partnerdeath)(1-error)(1-partnererror);
(* Singles due to switch, death, or error *)
single=Flatten[Position[Issingle,1]];
(*Selection & mutation process*)
Ndeath=Total[death];
Pdeath=Flatten[Position[death,1]];
newborn=RandomChoice[fitness->Range[n],Ndeath];
(* Newborn survival depends on the fitness of parents *)
newtraits=Table[{trait1[[newborn[[i]]]],trait2[[newborn[[i]]]]},{i,Ndeath}];
(* Newborn traits are inherited from parents *)
mutation=RandomInteger[BernoulliDistribution[mu],{Ndeath,2}];
(* Independent mutations occur on the 2 traits *)
newtraits=mutation*(3-2 newtraits)+newtraits;
If[Ndeath >0, 
trait1[[Pdeath]]=Transpose[newtraits][[1]];
trait2[[Pdeath]]=Transpose[newtraits][[2]],
trait1=trait1;
trait2=trait2];
type=trait1+trait2*2-2;
xt={N[Count[type,1]/n],N[Count[type,2]/n],N[Count[type,3]/n],N[Count[type,4]/n]};
xtlist=Append[xtlist,xt],
{step}];
xlist=Append[xlist,xtlist],
{run}]
];
xlist]
End[]

EndPackage[]
