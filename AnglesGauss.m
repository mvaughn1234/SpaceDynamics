% Obsv : [RA1, DEC1; RA2, DEC2; RA3, DEC3];
% R : [Site1,Site2,Site3];
% T : [T1; T2; T3]; with Ti = ti - t2
function r_set = AnglesGauss(obsv,R,T,mu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
RA = obsv(:,1)*pi/180;
DEC = obsv(:,2)*pi/180;
R1 = R(:,1); R2 = R(:,2); R3 = R(:,3);
T1 = T(1); T2 = T(2); T3 = T(3);

p1_unit = [cos(DEC(1))*cos(RA(1));...
    cos(DEC(1))*sin(RA(1));...
    sin(DEC(1))];

p2_unit = [cos(DEC(2))*cos(RA(2));...
    cos(DEC(2))*sin(RA(2));...
    sin(DEC(2))];

p3_unit = [cos(DEC(3))*cos(RA(3));...
    cos(DEC(3))*sin(RA(3));...
    sin(DEC(3))];

a1 = T3/(T3-T1); a3 = -T1/(T3-T1);
a1u = T3*((T3-T1)^2-T3^2)/(6*(T3-T1));
a3u = -T1*((T3-T1)^2-T1^2)/(6*(T3-T1));

L = [p1_unit,p2_unit,p3_unit];

M = inv(L)*R;

A = M(2,1)*a1-M(2,2)+M(2,3)*a3;
B = M(2,1)*a1u+M(2,3)*a3u;

E = dot(p2_unit,R2);

syms r2;
a = -(A^2+2*A*E+dot(R2,R2)); b = -2*mu*B*(A+E); c = -mu^2*B^2;
r2_set = solve(r2^8+a*r2^6+b*r2^3+c==0,r2);
r2 = vpa(r2_set(2));

u = mu/(r2^3);
p2 = A+B*u;

c1 = a1+a1u*u; c3 = a3+a3u*u;
c = [c1;-1;c3];
cp = -M*c;

p1 = cp(1)/c1;
p3 = cp(3)/c3;

r1 = p1*p1_unit + R1;
r2 = p2*p2_unit + R2;
r3 = p3*p3_unit + R3;
r_set = [r1,r2,r3];
end

