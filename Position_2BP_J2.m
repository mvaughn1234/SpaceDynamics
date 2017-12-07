function [t2,RV2] = Position_2BP_J2(initial_conds,t_range,tolerance,Re,J2,u)
    y_0 = [reshape(initial_conds,[6 1])];
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [t2,RV2] = ode45(@propagate_2BP_J2,t_range,y_0,options,Re,J2,u);
end

function out = propagate_2BP_J2(t,RV,Re,J2,u)
    out = [RV(4:6); -u*RV(1:3)/norm(RV(1:3))^3];
    out = out + [0;0;0;((-3*u*(Re^2)*J2)/(2*norm(RV(1:3))^5))*(1-5*(RV(3)^2)/norm(RV(1:3))^2)*RV(1);0;0];
    out = out + [0;0;0;0;((-3*u*(Re^2)*J2)/(2*norm(RV(1:3))^5))*(1-5*(RV(3)^2)/norm(RV(1:3))^2)*RV(2);0];
    out = out + [0;0;0;0;0;((-3*u*(Re^2)*J2)/(2*norm(RV(1:3))^5))*(3-5*(RV(3)^2)/norm(RV(1:3))^2)*RV(3)];
end