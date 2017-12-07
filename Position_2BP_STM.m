function [t2,RV2,STM] = Position_2BP_STM(initial_conds,t_range,tolerance,Re,J2,u)
    y_0 = [reshape(initial_conds,[6 1]); reshape(eye(6),[36 1])];
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [t2,out] = ode45(@propagate_2BP,t_range,y_0,options,Re,J2,u);
    RV2 = out(:,1:6);
    STM = out(:,7:6+36);
end

function dy = propagate_2BP(t,RV,Re,J2,u)
    dy = zeros(6+36,1);
    % For 2-Body Dynamics
    dy(1:3) = RV(4:6);
    dy(4:6) = (-u/norm(RV(1:3))^3)*RV(1:3);
    % For J2:
    dy(4) = dy(4) + ((-3*u*(Re^2)*J2)/(2*norm(RV(1:3))^5))*(1-5*(RV(3)^2)/norm(RV(1:3))^2)*RV(1);
    dy(5) = dy(5) + ((-3*u*(Re^2)*J2)/(2*norm(RV(1:3))^5))*(1-5*(RV(3)^2)/norm(RV(1:3))^2)*RV(2);
    dy(6) = dy(6) + ((-3*u*(Re^2)*J2)/(2*norm(RV(1:3))^5))*(3-5*(RV(3)^2)/norm(RV(1:3))^2)*RV(3);
    % For STM:
    a = -u/norm(RV(1:3))^3; 
    b = 3*u/norm(RV(1:3))^5; 
    Phi = reshape(RV(7:6+36,1),[6 6]); 
    A = [0, 0, 0, 1, 0, 0;
        0, 0, 0, 0, 1, 0;
        0, 0, 0, 0, 0, 1;
        a + b*RV(1)^2,b*RV(1)*RV(2),b*RV(1)*RV(3),0,0,0;
        b*RV(1)*RV(2),a + b*RV(2)^2,b*RV(2)*RV(3),0,0,0;
        b*RV(1)*RV(3),b*RV(2)*RV(3),a + b*RV(3)^2,0,0,0]*Phi;
    dy(7:6+36,1) = reshape(A,[36 1]); 
end