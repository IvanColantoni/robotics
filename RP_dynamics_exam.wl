(* ::Package:: *)

Clear
ClearAll
ClearAttributes
Q={q1[t],q2[t]};
dQ=D[Q,t];

p={(q2[t]-d)*Cos[q1[t]],(q2[t]-d)*Sin[q1[t]]};
J=D[p,{Q}];
v=D[p,t];
U=D[(q2[t]-d)*Sin[q1[t]],{Q}]
T1=1/2*I1*q1'[t]^2;
T2=1/2*v.v*m+1/2*I2*q1'[t]^2
T=T1+T2;
T=Expand[T];
Coll0=TrigReduce[T];
Coll1=Collect[Coll0,q1'[t]^2];
Coll2=Collect[Coll1,q2'[t]^2];
Coll3=Collect[Coll1,q2'[t]*q1'[t]]
M={{I1+I2+d^2 m-2 d m q2[t]+m q2[t]^2,0},{0,m}};
M.dQ;
dM=D[M,t];
MatrixForm[dM]
C1=(1/2)*(Table[D[M[[All,1]],Q[[i]]],{i,1,2}]+Transpose[Table[D[M[[All,1]],Q[[i]]],{i,1,2}]]-D[M,q1[t]]);
C2=(1/2)*(Table[D[M[[All,2]],Q[[i]]],{i,1,2}]+Transpose[Table[D[M[[All,2]],Q[[i]]],{i,1,2}]]-D[M,q2[t]]);

MatrixForm[C1]
MatrixForm[C2]

S1=C1.dQ;
S2=C2.dQ;
S={S1,S2};
MatrixForm[S]
Cq={dQ.S1,dQ.S2};
MatrixForm[Simplify[Cq]]
MatrixForm[Simplify[dM-2*S]]
St=Transpose[S];
MatrixForm[St.dQ]



