(* ::Package:: *)

(* symbolic computation of centrifugal and coriolis terms *)
M=  {{m1+m2+m3 ,-m2*d*Sin[q2]-m3*q3*Sin[q2], m3*Cos[q2]},
	{-m2*d*Sin[q2]-m3*q3*Sin[q2],m2*d^2+m3*q3^2+I2+I3, 0},
	{m3*Cos[q2],0 , m3}};
MatrixForm[M]
dM = {{0,-dq2*m2*d*Cos[q2]-m3*dq3*Sin[q2]-m3*q3*dq2*Cos[q2],-m3*dq2*Sin[q2]},
	{-dq2*m2*d*Cos[q2]-m3*dq3*Sin[q2]-m3*q3*dq2*Cos[q2],2*m3*q3*dq3,0},
	{-m3*dq2*Sin[q2],0,0}};
MatrixForm[dM]
Q ={q1,q2,q3};
dQ ={dq1,dq2,dq3};

C1=(1/2)*(Table[D[M[[All,1]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,1]],Q[[i]]],{i,1,3}]]-D[M,q1])
C2=(1/2)*(Table[D[M[[All,2]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,2]],Q[[i]]],{i,1,3}]]-D[M,q2])
C3=(1/2)*(Table[D[M[[All,3]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,3]],Q[[i]]],{i,1,3}]]-D[M,q3])
C11=dQ.C1.dQ;
C22=dQ.C2.dQ;
C33=dQ.C3.dQ;
Cq={C11,C22,C33};
MatrixForm[Cq]
C1S=(1/2)*Transpose[(Table[D[M[[All,1]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,1]],Q[[i]]],{i,1,3}]])].dQ;
(*C2S=(1/2)*(Table[D[M[[All,2]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,2]],Q[[i]]],{i,1,3}]]).dQ;*)
C2S=(1/2)*Transpose[(Table[D[M[[All,2]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,2]],Q[[i]]],{i,1,3}]])].dQ;
C3S=(1/2)*Transpose[(Table[D[M[[All,3]],Q[[i]]],{i,1,3}]+Transpose[Table[D[M[[All,3]],Q[[i]]],{i,1,3}]])].dQ;
CS={C1S,C2S,C3S};
MatrixForm[CS];
(*factorization of c (not unique), we can do it to obtain symmetric properties as: dM-2Cfac symmetric*)

Cfac = {dQ.C1,dQ.C2,dQ.C3};
V1= dM-2*Cfac;
(*V2 = dM-2*CS;*)
MatrixForm[Cfac]

MatrixForm[Simplify[V1]]
MatrixForm[Simplify[V2]]
(*potential energy of the system U, Gravitational term G, Hessian of Potential energy G2*)
U1=0;
U2=m2*g0*d*Sin[q2];
U3=m3*g0*q3*Sin[q2];
U= U1+U2+U3;
G=Table[D[U,Q[[i]]],{i,1,3}];
MatrixForm[G]

G2=Table[D[G,Q[[i]]],{i,1,3}];
MatrixForm[G2.G2]








