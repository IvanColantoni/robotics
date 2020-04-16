(* ::Package:: *)

Q={q1[t],q2[t]};
dQ=D[Q,t];
M={{a1+2*a2*Cos[q2[t]],a3+a2*Cos[q2[t]]},{a3+a2*Cos[q2[t]],a3}};
dM=D[M,t];
MatrixForm[dM]
C1=(1/2)*(Table[D[M[[All,1]],Q[[i]]],{i,1,2}]+Transpose[Table[D[M[[All,1]],Q[[i]]],{i,1,2}]]-D[M,q1[t]]);
C2=(1/2)*(Table[D[M[[All,2]],Q[[i]]],{i,1,2}]+Transpose[Table[D[M[[All,2]],Q[[i]]],{i,1,2}]]-D[M,q2[t]]);
C11=dQ.C1.dQ;
C22=dQ.C2.dQ;
Sa1=C1.dQ;
Sa2=C2.dQ;
S1={Sa1,Sa2};
MatrixForm[S1]
Cq={C11,C22};
MatrixForm[Simplify[Cq]]
MatrixForm[Simplify[S1.dQ]]
MatrixForm[Simplify[dM-2*S1]]
S2={{-2*a2*q2'[t]*Sin[q2[t]],-a2*q2'[t]*Sin[q2[t]]},{a2*q1'[t]*Sin[q2[t]],0}};
MatrixForm[S2]
MatrixForm[dM-2*S2]



