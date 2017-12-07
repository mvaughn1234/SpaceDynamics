function [tesc, tesc_true,t2,RV2] = low_thrust_escape_time(initial_conds,thrust,u,tol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
r0 = norm(initial_conds(1:3));
v0 = norm(initial_conds(4:6));
aT = thrust;
tesc = (v0/aT)*(1-((20*(aT^2)*(r0^2))/(v0^4))^(1/8));
[t2 RV2] = Position_2BP(initial_conds,aT,[0 10*tesc],tol,u);

length = size(RV2(:,1));
length = length(1);
tesc_true = tesc;
first = true; first_t = true;

for i = [1:1:length]
    rcurr = norm(RV2(i,1:3));
    vcurr = norm(RV2(i,4:6));
    E = (vcurr^2)/2 - u/rcurr;
    if(E >= 0 && first)
        first = false;
        tesc_true = t2(i);
        Esc = E;
        resc = rcurr;
        vesc = vcurr;
    end
    if(t2(i) >= tesc && first_t)
        first_t = false;
        rcurr = norm(RV2(i,1:3));
        vcurr = norm(RV2(i,4:6));
        E_est = (vcurr^2)/2 - u/rcurr;
    end
end
end

