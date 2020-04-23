%mobility analysis
syms q1 q2 q3;
q=[q1 q2 q3].';
p=[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
    sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
Jac=jacobian(p,q);

fprintf('direct kinematics of 3-R planar robot with unitary links\n');
fprintf('Jac=\n');
%evaluating the Jacobian in that configuration
J=(subs(Jac,{q1,q2,q3},{0,pi/2,pi/4}));
disp(vpa(J,5));
pause;

%I'm asking if all the directions of the end-effector are feasible from
%that configuration applying any joint velocity. I have to know the rank of the matrix J computed as
%above. This gives me the dimension of the subspace of velocities
%achievable from that configuration.
rankJ=rank(J);
fprintf('rankJ:\n');
disp(rankJ);
pause;

%Which is a base for that subspace? range(J) gives me a set of unitary norm
%vector that, that set is a base. 
rangeJ=range(J);
fprintf('rangeJ:\n');
disp(rangeJ);
pause;

%I'm asking if there are velocities of the joints that allow the
%end-effector to be fixed in the task (cartesian) space. The command
%null(J) compute a base for the null space of J. This means that the values
%of the (unitary norm) vector kerJ are the joint velocties which make the
%end-effector to have zero velocity. rankj+dim(kerJ)=m
kerJ=null(J);
fprintf('kerJ:\n')
disp(vpa(kerJ,4));
pause;

%dual problem with Jacobian transpose
JT=J.';
fprintf('JacobianTranspose:\n');
disp(vpa(JT,4));

%rankJ=rankJT

%null-space of JT: I'm checking if there are possible directions/values of 
%applicable forces from the task space which don't need to be balanced with
%torques/forces in the joint space. But hence the matrix J has full rank
%(rank=m) his transpose JT won't generate a null space different from the zero vector {0}.
% Anyway this could be checked with the null space of JT.
kerJT=null(JT);
fprintf('kerJTranspose:\n')
disp(vpa(kerJT,4));
pause;

%check on the dimensions of supspaces:
fprintf('dim(RangeJ)+dim(KerJT)= n ?\n');
disp(rank(JT)+rank(kerJT));

pause
fprintf('dim(RangeJT)+dim(KerJ)= m?\n');
disp(rank(JT)+rank(kerJ));

