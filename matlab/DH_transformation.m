function[T]=DH_transformation(a,alpha,d,teta)
%angoli in radianti
T = [cos(teta) sin(teta)*(-cos(alpha)) sin(alpha)*sin(teta) a*cos(teta);
    sin(teta) cos(alpha)*cos(teta) -sin(alpha)*cos(teta) a*sin(teta);
    0 sin(alpha) cos(alpha) d;
    0 0 0 1];
 disp(T) 
end 