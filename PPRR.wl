(* ::Package:: *)

Q={q1[t],q2[t],q3[t],q4[t]};
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
Simplify[T];
Collect[T,q1'[t]^2]
Simplify[Collect[T,q2'[t]^2]]
Collect[T,q3'[t]^2]
Collect[T,q4'[t]^2]




