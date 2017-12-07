%input = [a;e;i;OMEGA;omega;nu]

clear all; close all; clc; format long;

u = 3.986*10^5;
input = [40000;0.2;140;0;90;180];

elementsToPQW(input,u);
new_input = elementsToECI(input,u).';
period = 2*pi*sqrt((input(1)^3)/(u));
trange = [1:period];

[t,RV] = Position_2BP(new_input,trange,10^-13,u);
initial_spec_e = (norm(RV(1,4:6))^2)/2-(u/norm(RV(1,1:3)));

for i = trange
    spec_e_Dev(i) = (norm(RV(i,4:6))^2)/2-(u/norm(RV(i,1:3))) - initial_spec_e;
end

% plot3(RV(:,1),RV(:,2),RV(:,3))
% axis equal
% 
% figure
% subplot(3,1,1)
% plot(t,RV(:,1))
% 
% subplot(3,1,2)
% plot(t,RV(:,2))
% 
% subplot(3,1,3)
% plot(t,RV(:,3))
% 
% figure
% plot(t,spec_e_Dev)
