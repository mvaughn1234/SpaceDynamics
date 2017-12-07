clear all; close all; clc;

u = 3.986*(10^5);
tolerance = 10^-13;
trange = [1:30000];
input = [5203.12878457022 2539.18526782417 4387.98384076804 -5.73055171828814 1.23647597198147 6.07959326945700];

[t,RV] = Position_2BP(input,trange,tolerance,u);

for i = trange
   eles(i,1:6) = orbital_elements(RV(i,1:3).',RV(i,4:6).',u);
end
% 
% subplot(2,3,1);
% plot(t,eles(:,1))
% 
% subplot(2,3,2);
% plot(t,eles(:,2))
% 
% subplot(2,3,3);
% plot(t,eles(:,3))
% 
% subplot(2,3,4);
% plot(t,eles(:,4))
% 
% subplot(2,3,5);
% plot(t,eles(:,5))
% 
% subplot(2,3,6);
% plot(t,eles(:,6))

function set = orbital_elements(r,v,u)

h = cross(r,v);
n = cross([0;0;1],h);
e = (1/u)*((norm(v)^2-u/norm(r))*r-(dot(r,v))*v);
p = (norm(h)^2)/u;
a = p/(1-norm(e)^2);
OMEGA = acos(dot(n,[1;0;0])/norm(n))*180/pi;
omega = acos(dot(n,e)/(norm(n)*norm(e)))*180/pi;
i = acos(dot(h,[0;0;1])/norm(h))*180/pi;
nu = acos(dot(r,e)/(norm(r)*norm(e)))*180/pi;

set = [a,norm(e),i,OMEGA,omega,nu]; 
end