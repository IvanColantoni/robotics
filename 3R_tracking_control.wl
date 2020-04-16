(* ::Package:: *)

Q={q1[t],q2[t],q3[t]};
p={Cos[q1[t]]+Cos[q2[t]]+Cos[q3[t]],Sin[q1[t]]+Sin[q2[t]]+Sin[q3[t]]};
pd={1+2*Sin[3*t],2+Cos[3*t+Pi/2]};
pd1=D[pd,t];
pd2=D[pd1,t];
pd3=D[pd2,t];

J=D[p,{Q}];
J1=D[J,t];
J2=D[J1,t];
(*Jn={{0,-1,-1},{1, 0, 0}};
Jps= PseudoInverse[Jn];
Jps.{6,-3};
J1n={{3,0,0},{0,3,3}};
J2n={{18,9,9},{9,-4.5^2,-4.5^2}}
Jps.({-54,27}-2*{{3,0,0},{0,3,3}}.{-18,4.5,4.5}-J2n.{-3,-3,-3})
Jps.J1n.{-3,-3,-3}
MatrixForm[J];
MatrixForm[J1]
MatrixForm[J2]
pd1
pd2
pd3*)

Jn={{-1,0,0},{0,1,1}};
Jn.



