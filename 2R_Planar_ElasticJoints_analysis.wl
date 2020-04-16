(* ::Package:: *)

Q={q1[t],q2[t],q3[t],q4[t]};
dQ={q1'[t],q2'[t],q3'[t],q4'[t]};
dQ2={q1''[t],q2''[t],q3''[t],q4''[t]};

pc1={d1*Cos[q1[t]],d1*Sin[q1[t]],0};
pc2={l1*Cos[q1[t]]+d2*Cos[q2[t]+q1[t]],l1*Sin[q1[t]]+d2*Sin[q1[t]+q2[t]],0};
v1=D[pc1,t];
v2=D[pc2,t];

T1=(1/2)*m1*v1.v1+(1/2)*I1*q1'[t]^2;
T2=(1/2)*m2*v2.v2+(1/2)*I2*(q1'[t]+q2'[t])^2;
T1m=(1/2)*I1m*q3'[t]^2;
T2m=(1/2)*I2m*(q1'[t]+q4'[t])^2+(1/2)*m2m*v1.v1;
T=T1+T2+T1m+T2m;
T=Expand[T];
Coll0=TrigReduce[T];
Coll1=Collect[Coll0,q1'[t]^2];
Coll2=Collect[Coll1,q2'[t]^2];
Coll3=Collect[Coll2,q3'[t]^2];
Coll4=Collect[Coll3,q4'[t]^2];
Coll5=Collect[Coll4,q1'[t]*q2'[t]];
Coll6=Collect[Coll5,q1'[t]*q3'[t]];
Coll7=Collect[Coll6,q1'[t]*q4'[t]];
Coll8=Collect[Coll7,q2'[t]*q3'[t]];
Coll9=Collect[Coll8,q2'[t]*q4'[t]];
Coll10=Collect[Coll9,q3'[t]*q4'[t]];
M={{a3+a5+a2+2*a4*Cos[q2[t]],a3+a4*Cos[q2[t]],0,a2},{a3+a4*Cos[q2[t]],a3,0,0},{0,0,a1,0},{a2,0,0,a2}};
MatrixForm[M]
MOD=M.dQ2
A={a1,a2,a3,a4,a5};
Y=CoefficientArrays[MOD,A]
Y[[2]]
TableForm[Simplify[Normal[%]]]





