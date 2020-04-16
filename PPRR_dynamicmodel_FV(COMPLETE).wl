(* ::Package:: *)

Q={q1[t],q2[t],q3[t],q4[t]};
dQ={q1'[t],q2'[t],q3'[t],q4'[t]};
dQ2={q1''[t],q2''[t],q3''[t],q4''[t]};
(*velocity of the centre of mass of each link*)
p1={q1[t]-d1,0,0};
p2={q1[t],q2[t]-d2,0};
p3={q1[t]+d3*Cos[q3[t]],q2[t]+d3*Sin[q3[t]],0};
p4={q1[t]+l3*Cos[q3[t]]+d4*Cos[q3[t]+q4[t]],q2[t]+l3*Sin[q3[t]]+d4*Sin[q3[t]+q4[t]],0};
v1=D[p1,t]
v2=D[p2,t]
v3=D[p3,t]
v4=D[p4,t]

T1=(1/2)*m1*v1.v1;
T2=(1/2)*m2*v2.v2;
T3=(1/2)*m3*v3.v3+(1/2)*I3*q3'[t]^2;
T4=(1/2)*m4*v4.v4+(1/2)*I4*(q3'[t]+q4'[t])^2;
T= T1+T2+T3+T4;
T=Expand[T];
(*CoefficientArrays[T,{q1'[t]^2,q2'[t]^2,q3'[t]^2,q4'[t]^2,q1'[t]*q2'[t],q1'[t]*q3'[t]}]*)
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
Coll10=Collect[Coll9,q3'[t]*q4'[t]]



M={{a1,0,-a6*Sin[q3[t]]-a5*Sin[q3[t]+q4[t]],-a5*Sin[q3[t]+q4[t]]},
{0,a2,a6*Cos[q3[t]]+a5*Cos[q3[t]+q4[t]],a5*Cos[q3[t]+q4[t]]},
{-a6*Sin[q3[t]]-a5*Sin[q3[t]+q4[t]],a6*Cos[q3[t]]+a5*Cos[q3[t]+q4[t]],a3+2*a5*l3*Cos[q4[t]],a4+a5*l3*Cos[q4[t]]},
{-a5*Sin[q3[t]+q4[t]],a5*Cos[q3[t]+q4[t]],a4+a5*l3*Cos[q4[t]],a4}};
MatrixForm[M]
C1=(1/2)*(Table[D[M[[All,1]],Q[[i]]],{i,1,4}]+Transpose[Table[D[M[[All,1]],Q[[i]]],{i,1,4}]]-D[M,q1[t]]);
C2=(1/2)*(Table[D[M[[All,2]],Q[[i]]],{i,1,4}]+Transpose[Table[D[M[[All,2]],Q[[i]]],{i,1,4}]]-D[M,q2[t]]);
C3=(1/2)*(Table[D[M[[All,3]],Q[[i]]],{i,1,4}]+Transpose[Table[D[M[[All,3]],Q[[i]]],{i,1,4}]]-D[M,q3[t]]);
C4=(1/2)*(Table[D[M[[All,4]],Q[[i]]],{i,1,4}]+Transpose[Table[D[M[[All,4]],Q[[i]]],{i,1,4}]]-D[M,q4[t]]);
MatrixForm[C1];
MatrixForm[C2];
MatrixForm[C3];
MatrixForm[C4];
C11=dQ.C1.dQ;
C22=dQ.C2.dQ;
C33=dQ.C3.dQ;
C44=dQ.C4.dQ;
Cq={C11,C22,C33,C44};
MatrixForm[Simplify[Cq]]
(*Potential energy*)
U1= 0;
U2=(q2[t]-d2)*g0*m2;
U3=(q2[t]+d3*Sin[q3[t]])*g0*m3;
U4=(q2[t]+l3*Sin[q3[t]]+d4*Sin[q3[t]+q4[t]])*g0*m4;
U=U1+U2+U3+U4;
(*vector G(q)
G=Table[D[U,Q[[i]]],{i,1,4}];
MatrixForm[G];*)
a2=.
a={a1,a2,a3,a4,a5,a6,a7,a8,a9,a10};
G={0,g0*a2,a6*g0*Cos[q3[t]]+a5*g0*Cos[q3[t]+q4[t]],a5*g0*Cos[q3[t]+q4[t]]};
Fv={a7*q1'[t],a8*q2'[t],a9*q3'[t],a10*q4'[t]};
u={u1,u2,u3,u4};

u=M.dQ2+Cq+G+Fv;
Y=CoefficientArrays[u,a]
Y[[2]]
TableForm[Simplify[Normal[%]]]




