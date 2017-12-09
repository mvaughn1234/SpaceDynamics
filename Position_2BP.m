function [t2,RV2] = Position_2BP(initial_conds,t_range,tolerance,u)
    y_0 = [reshape(initial_conds,[6 1])];
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [t2,RV2] = ode45(@propagate_2BP,t_range,initial_conds,options,u);
end

function out = propagate_2BP(t,RV,u)
    out = [RV(4:6); -u*RV(1:3)/norm(RV(1:3))^3];
end