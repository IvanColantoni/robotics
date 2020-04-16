(* ::Package:: *)

Q={q1[t],q2[t],q3[t],q4[t]};
Pe={l*Cos[q1[t]]+l*Cos[q1[t]+q2[t]]+l*Cos[q1[t]+q2[t]+q3[t]]+l*Cos[q1[t]+q2[t]+q3[t]+q4[t]],
	l*Sin[q1[t]]+l*Sin[q1[t]+q2[t]]+l*Sin[q1[t]+q2[t]+q3[t]]+l*Sin[q1[t]+q2[t]+q3[t]+q4[t]],
	0};
J=Table[D[Pe[[i]],{Q}],{i,1,3}]
TrigReduce[J];
MatrixForm[%]



