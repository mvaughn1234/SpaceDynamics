clear all; close all; clc;
u = 3.986E5;
tol = 1E-13;
i = 1;

r = [8000 0 0];
v0 = sqrt(u/8000);
v = [0 v0 0];
i_c = [r v];

thrust = linspace(1E-5,1E-3,100);
for aT = thrust
    tesc(i,1) = aT;
    tesc_true(i,1) = aT;
    [tesc(i,2), tesc_true(i,2)] = low_thrust_escape_time(i_c,aT,u,tol);
    i = i+1;
end
figure;
hold on;
plot(tesc(:,1),tesc(:,2),'--black');
plot(tesc_true(:,1),tesc_true(:,2),'black');
xlabel('Specific Thrust [kN/kg]');
ylabel('Time to Escape [sec] (log)');
legend('t_e_s_c','t_e_s_c_\__t_r_u_e','Location','NorthEast');
set(gca,'yscale','log');
grid on;