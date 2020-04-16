(* ::Package:: *)

Clear; 
ClearAll;
ClearAttributes;
Q={q1[t],q2[t],q3[t]};
dQ=D[Q,t];

p1={q1[t],0,0};
p2={q1[t]+d*Cos[q2[t]],d*Sin[q2[t]],0};
p3={q1[t]+q3[t]*Cos[q2[t]],q3[t]*Sin[q2[t]],0};

v1=D[p1,t];
v2=D[p2,t];
v3=D[p3,t];

T1=1/2*m1*v1.v1;
T2=1/2*m2*v2.v2+ 1/2*I2*q2'[t]^2;
T3=1/2*m2*v3.v3+ 1/2*I3*q2'[t]^2;
T=T1+T2+T3;
T=Expand[T];
Coll0=TrigReduce[T];
Coll1=Collect[Coll0,q1'[t]^2];
Coll2=Collect[Coll1,q2'[t]^2];
Coll3=Collect[Coll2,q3'[t]^2];
Coll4=Collect[Coll3,q1'[t]*q2'[t]];
Coll5=Collect[Coll4,q1'[t]*q3'[t]];
Coll6=Collect[Coll5,q2'[t]*q3'[t]]

M={{m1+m2,-Sin[q2[t]]*(d*m2+m2*q3[t]),m2*Cos[q2[t]]},
	{-Sin[q2[t]]*(d*m2+m2*q3[t]),I2+I3+d^2*m2+m2*q3[t]^2,0},
	{m2*Cos[q2[t]],0,m2}};

Mp={{a1+a2,-(a4+a2*q3[t])*Sin[q2[t]],a2*Cos[q2[t]]},
	{-(a4+a2*q3[t])*Sin[q2[t]],a3+a5+a2*q3[t]^2,0},
	{a2*Cos[q2[t]],0,a2}};
	
dM=D[Mp,t];
MatrixForm[dM]
C1=(1/2)*(Table[D[Mp[[All,1]],Q[[i]]],{i,1,3}]+Transpose[Table[D[Mp[[All,1]],Q[[i]]],{i,1,3}]]-D[Mp,q1[t]]);
C2=(1/2)*(Table[D[Mp[[All,2]],Q[[i]]],{i,1,3}]+Transpose[Table[D[Mp[[All,2]],Q[[i]]],{i,1,3}]]-D[Mp,q2[t]]);
C3=(1/2)*(Table[D[Mp[[All,3]],Q[[i]]],{i,1,3}]+Transpose[Table[D[Mp[[All,3]],Q[[i]]],{i,1,3}]]-D[Mp,q3[t]]);

MatrixForm[C1];
MatrixForm[C2];
MatrixForm[C3];

S1=C1.dQ;
S2=C2.dQ;
S3=C3.dQ;

S={S1,S2,S3};
Cq={dQ.S1,dQ.S2,dQ.S3}

MatrixForm[S]
MatrixForm[Simplify[Cq]]

MatrixForm[dM-2*S];
Simplify[dM-2*S]

U=m2*g0*d*Sin[q2[t]]+m3*g0*q3[t]*Sin[q2[t]];
G=Table[D[U,Q[[i]]],{i,1,3}]
H=D[G,{Q}]


MatrixForm[H.H]

Eigenvalues[H.H]
Simplify[TrigReduce[%]]



