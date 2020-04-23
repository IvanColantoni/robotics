%This program's aim is to find vci and wci relatade to each joint of the
%manipulator. Moving frame algorithm is used to achieve that scope. First,
%in the case we would use DH-transformation matrices.. (then) 
%For the moment, given related properties of each joint (revolute or
%prismatic) we assign the S[0,..., i-1]  vector which contains the values
%(0 or 1) in the case that the joint is revolute or prismatic respectively.
%

%function[T]=DH_transformation(a,alpha,d,teta)
%angoli in radianti
%T = [cos(teta) sin(teta)*(-cos(alpha)) sin(alpha)*sin(teta) a*cos(teta);
%    sin(teta) cos(alpha)*cos(teta) -sin(alpha)*cos(teta) a*sin(teta);
%    0 sin(alpha) cos(alpha) d;
%    0 0 0 1];
% disp(T) 
% end 

%number of joints/generalized variables
n = 2;
S = zeros(1,2,'sym');
W= sym('W',[3 3]);
W(:,1) = zeros(1,3,'sym');

%V= cell(n+1);
%linear velocity 2nd row
V= sym('V',[3 3]);
V(:,1) = zeros(1,3,'sym');
%sym('Q',[1 n]);

syms m1 m2 v1 v2 q1 q2 dq1 dq2 w1 w2 l1 l2 d1 d2 ;
L = sym('L', [3 2]);
L(:,1)= [l1*cos(q1), l1*cos(q1), 0];
L(:,2) =[l2*cos(q2), l2*cos(q2), 0];
%Q(1)= q1;
Q= sym('Q',[2 2]);
Q(1,1)= q1;
Q(2,1) = q2;
Q(2,2) = dq2;
Q(1,2) = dq1;

%dQ = cell(n);
%Q{2,1}= zeros(1,'sym');
%Q{2,2} = dq1;
%Q{2,3} = dq2;

R = cell(n+1);

R{1,1} = [ cos(q1), sin(q1), 0; -sin(q1), cos(q1),0;  0,       0, 1];
R{1,2} = [cos(q2), sin(q2),  0;
      -sin(q1), cos(q1), 0;
             0,       0, 1];
         
for i=1:n
    W(:,i+1) = transpose(R{1,i})*(W(:,i) + (1-S(i))*Q(i,2)*sym([0,0,1],'f')');
    V(:,i+1) = transpose(R{1,i})*((V(:,i) +S(i)*(Q(i,2)*sym([0,0,1],'f')')+cross(W(:,i+1),L(:,i))));
end


disp(W);
disp(V);


