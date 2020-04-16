(* ::Package:: *)

(*This wolphram code provides the symbolic expression of a RRR planar manipulator with respect to generalized variables Q={q1[t],q2[t],q3[t]},
and their relative derivatives wrt parameter [t]. *)


ClearAll
Q={q1[t],q2[t],q3[t]};
(*velocity of the centre of mass of each link*)
P={Cos[q1[t]]+Cos[q1[t]+q2[t]]+Cos[q1[t]+q2[t]+q3[t]],Sin[q1[t]]+Sin[q1[t]+q2[t]]+Sin[q1[t]+q2[t]+q3[t]],0};
J=D[P,{Q}];
MatrixForm[J]
v1={-0.25*q1'[t]*Sin[q1[t]],-0.25*q1'[t]*Cos[q1[t]],0};
v2={-(0.5*q1'[t])*Sin[q1[t]]-0.25*q2'[t]*Sin[q2[t]],(0.5*q1'[t])*Cos[q1[t]]+0.25*q2'[t]*Cos[q2[t]],0};
v3={-(0.5*q1'[t])*Sin[q1[t]]-(0.5*q2'[t])*Sin[q2[t]]-0.25*q3'[t]*Sin[q3[t]],(0.5*q1'[t])*Cos[q1[t]]+(0.5*q2'[t])*Cos[q2[t]]+0.25*q3'[t]*Cos[q3[t]],0};(*dQ={dq1,q2'[t],q3'[t]};*)
dQ=D[Q,t];
T1= 0.5*30*v1.v1+0.5*(1/12)*30*(0.5)^2*dq1^2;
T2= 0.5*20*v2.v2+0.5*(1/12)*20*(0.5)^2*(dq1+q2'[t])^2;
T3= 0.5*10*v3.v3+0.5*(1/12)*10*(0.5)^2*(dq1+q2'[t]+q3'[t])^2;
T= T1 + T2 +T3;
T = Simplify[T1+T2+T3];
M={{2*5.3125, 5*Cos[q1[t]-q2[t]]*0.625,0.20833+1.25*Cos[q1[t]-q3[t]]},
{5*Cos[q1[t]-q2[t]]*0.625, 2.1875*2 , 0.20833+1.25*Cos[q2[t]-q3[t]]},
{0.20833+1.25*Cos[q1[t]-q3[t]],0.20833+1.25*Cos[q2[t]-q3[t]],0.4166666666666667`}};
MatrixForm[M];
dM=D[M,t];
MatrixForm[dM]
C1=(1/2)*(Table[D[M[[All,1]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,1]],Q[[i]]],{i,1,3}]]-D[M,q1[t]]);
C2=(1/2)*(Table[D[M[[All,2]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,2]],Q[[i]]],{i,1,3}]]-D[M,q2[t]]);
C3=(1/2)*(Table[D[M[[All,3]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,3]],Q[[i]]],{i,1,3}]]-D[M,q3[t]]);
(*C11=dQ.C1.dQ;
C22=dQ.C2.dQ;
C33=dQ.C3.dQ;
Cq={C11,C22,C33};
MatrixForm[Cq];*)
Cfac = {dQ.C1,dQ.C2,dQ.C3};
MatrixForm[Cfac]
Dimensions[Cfac]
V=dM-2*Cfac
Dimensions[dM]
Dimensions[V]
MatrixForm[Simplify[V]]

Dimensions[Cfac]
U1=30*9.81*0.25*Sin[q1[t]];
U2=20*9.81*(0.5*Sin[q1[t]]+0.25*Sin[q2[t]]);
U3=10*9.81*(0.5*Sin[q1[t]]+0.5*Sin[q2[t]]+0.25*Sin[q3[t]]);
U=U1+U2+U3;

G=Table[D[U,Q[[i]]],{i,1,3}];
MatrixForm[G]




(* ::InheritFromParent:: *)
(*0.20833333333333331` (dq1+q2'[t])^2+10.` ((0.5` dq1 Cos[q1]+0.25` q2'[t] Cos[q2])^2+(-0.5` dq1 Sin[q1]-0.25` q2'[t] Sin[q2])^2)*)


(* ::InheritFromParent:: *)
(*0.10416666666666666` (dq1+q2'[t]+q3'[t])^2+5.` ((0.5` dq1 Cos[q1]+0.5` q2'[t] Cos[q2]+0.25` q3'[t] Cos[q3])^2+(-0.5` dq1 Sin[q1]-0.5` q2'[t] Sin[q2]-0.25` q3'[t] Sin[q3])^2)*)


(* ::InheritFromParent:: *)
(*0.3125` dq1^2+0.20833333333333331` (dq1+q2'[t])^2+0.10416666666666666` (dq1+q2'[t]+q3'[t])^2+15.` (0.0625` dq1^2 Cos[q1]^2+0.0625` dq1^2 Sin[q1]^2)+10.` ((0.5` dq1 Cos[q1]+0.25` q2'[t] Cos[q2])^2+(-0.5` dq1 Sin[q1]-0.25` q2'[t] Sin[q2])^2)+5.` ((0.5` dq1 Cos[q1]+0.5` q2'[t] Cos[q2]+0.25` q3'[t] Cos[q3])^2+(-0.5` dq1 Sin[q1]-0.5` q2'[t] Sin[q2]-0.25` q3'[t] Sin[q3])^2)*)


(* ::InheritFromParent:: *)
(*5.3125` dq1^2+0.625` dq1 q2'[t]+2.1875` q2'[t]^2+0.20833333333333334` dq1 q3'[t]+0.20833333333333334` q2'[t] q3'[t]+0.4166666666666667` q3'[t]^2+5.` dq1 q2'[t] Cos[q1-q2]+1.25` dq1 q3'[t] Cos[q1-q3]+1.25` q2'[t] q3'[t] Cos[q2-q3]*)


(* ::InheritFromParent:: *)
(*\!\( *)
(*TagBox[*)
(*RowBox[{"(", "", *)
(*TagBox[GridBox[{*)
(*{*)
(*RowBox[{"0.`", " ", "+", *)
(*RowBox[{*)
(*RowBox[{*)
(*RowBox[{"q2", "'"}], "[", "t", "]"}], " ", *)
(*RowBox[{"(", *)
(*RowBox[{"0.`", " ", "+", *)
(*RowBox[{"3.125`", " ", *)
(*RowBox[{*)
(*RowBox[{"q2", "'"}], "[", "t", "]"}], " ", *)
(*RowBox[{"Sin", "[", *)
(*RowBox[{"q1", "-", "q2"}], "]"}]}]}], ")"}]}], "+", *)
(*RowBox[{*)
(*RowBox[{*)
(*RowBox[{"q3", "'"}], "[", "t", "]"}], " ", *)
(*RowBox[{"(", *)
(*RowBox[{"0.`", " ", "+", *)
(*RowBox[{"1.25`", " ", *)
(*RowBox[{*)
(*RowBox[{"q3", "'"}], "[", "t", "]"}], " ", *)
(*RowBox[{"Sin", "[", *)
(*RowBox[{"q1", "-", "q3"}], "]"}]}]}], ")"}]}]}]},*)
(*{*)
(*RowBox[{"0.`", " ", "+", *)
(*RowBox[{"dq1", " ", *)
(*RowBox[{"(", *)
(*RowBox[{"0.`", " ", "-", *)
(*RowBox[{"3.125`", " ", "dq1", " ", *)
(*RowBox[{"Sin", "[", *)
(*RowBox[{"q1", "-", "q2"}], "]"}]}]}], ")"}]}], "+", *)
(*RowBox[{*)
(*RowBox[{*)
(*RowBox[{"q3", "'"}], "[", "t", "]"}], " ", *)
(*RowBox[{"(", *)
(*RowBox[{"0.`", " ", "+", *)
(*RowBox[{"1.25`", " ", *)
(*RowBox[{*)
(*RowBox[{"q3", "'"}], "[", "t", "]"}], " ", *)
(*RowBox[{"Sin", "[", *)
(*RowBox[{"q2", "-", "q3"}], "]"}]}]}], ")"}]}]}]},*)
(*{*)
(*RowBox[{"0.`", " ", "+", *)
(*RowBox[{"dq1", " ", *)
(*RowBox[{"(", *)
(*RowBox[{"0.`", " ", "-", *)
(*RowBox[{"1.25`", " ", "dq1", " ", *)
(*RowBox[{"Sin", "[", *)
(*RowBox[{"q1", "-", "q3"}], "]"}]}]}], ")"}]}], "+", *)
(*RowBox[{*)
(*RowBox[{*)
(*RowBox[{"q2", "'"}], "[", "t", "]"}], " ", *)
(*RowBox[{"(", *)
(*RowBox[{"0.`", " ", "-", *)
(*RowBox[{"1.25`", " ", *)
(*RowBox[{*)
(*RowBox[{"q2", "'"}], "[", "t", "]"}], " ", *)
(*RowBox[{"Sin", "[", *)
(*RowBox[{"q2", "-", "q3"}], "]"}]}]}], ")"}]}]}]}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Column], "", ")"}],*)
(*Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)*)
