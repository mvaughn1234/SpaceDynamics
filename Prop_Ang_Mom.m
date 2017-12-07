function [t,wnew] = Prop_Ang_Mom(I,w,t,tolerance)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
[t,wnew] = ode45(@odeFun,t,w,options,I);
end

function wdot = odeFun(t,w,I)
wdot = [-(I(3,3)-I(2,2))*w(2,1)*w(3,1)/I(1,1);
    -(I(1,1)-I(3,3))*w(3,1)*w(1,1)/I(2,2);
    -(I(2,2)-I(1,1))*w(1,1)*w(2,1)/I(3,3);];
end