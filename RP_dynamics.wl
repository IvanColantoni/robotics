(* ::Package:: *)

Q={q1[t],q2[t]};
dQ=D[Q,t];
p1={d1*Cos[q1[t]],d1*Sin[q1[t]]};
p2={l1*Cos[q1[t]]-(q2[t]-d2)*Sin[q1[t]],l1*Sin[q1[t]]+(q2[t]-d2)*Cos[q1[t]]};

v1=D[p1,t];
v2=D[p2,t];

T1=1/2*v1.v1*m1+1/2*I1*q1'[t]^2;
T2=1/2*v2.v2*m2+1/2*I2*q1'[t]^2;

T=T1+T2;
T=Expand[T];
Coll0=TrigReduce[T];
Coll1=Collect[Coll0,q1'[t]^2];
Coll2=Collect[Coll1,q2'[t]^2];
Coll3=Collect[Coll1,q2'[t]*q1'[t]]


M={{I1+I2+d1^2*m1+d2^2*m2+l1^2*m2-2*d2*m2*q2[t]+m2*q2[t]^2,l1*m2},{l1*m2,m2}};
dM=D[M,t];
MatrixForm[M]
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

SS2={{-2*m2*(d2+q2[t])q2'[t],0},{m2*(d2+q2),0}};
MatrixForm[Simplify[SS2.dQ]]
MatrixForm[Simplify[dM-2*SS2]]



U=m1*g0*d1*Sin[q1[t]]+m2*(l1*Sin[q1[t]]+(q2[t]-d1)*Cos[q1[t]]);
G=Simplify[D[U,{Q}]]

P=MatrixForm[M.dQ]

JK1=D[p1,{Q}];
Eigensystem[JK1]

JK2=D[p2,{Q}]
JK2={{-l1 Sin[q1[t]],-Sin[q1[t]]},{l1 Cos[q1[t]],Cos[q1[t]]}}
Eigensystem[JK2]



