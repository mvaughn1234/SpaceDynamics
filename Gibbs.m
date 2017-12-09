% r_set : [r1,r2,r3]
function v1 = Gibbs(r_set,mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r1 = r_set(:,1); r1n = norm(r1);
r2 = r_set(:,2); r2n = norm(r2);
r3 = r_set(:,3); r3n = norm(r3);

D = cross(r1,r2) + cross(r2,r3) + cross(r3,r1);
N = r3n*cross(r1,r2) + r1n*cross(r2,r3) + r2n*cross(r3,r1);
S = (r2n-r3n)*r1 + (r3n-r1n)*r2 + (r1n-r2n)*r3;
Nn = norm(N); Dn = norm(D); Sn = norm(S);
v1 = sqrt(mu/(Nn*Dn))*(cross(D,(r1/r1n))+S);
end