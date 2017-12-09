%% ENAE 441 Final Project
% Mathew Vaughn
clear all; close all; clc;
format long; format compact
Re = 6371; % km
mu = 3.986E5;

obsv1 = [21,14,0,-28,0,0;...
    21,15,0,-15,0,0;...
    21,15,0,8,30,0;...
    21,15,0,47,15,0;...
    21,10,0,81,0,0;...
    9,30,0,81,0,0];

obsv2 = [9,16,35,82,20,35;...
    1,14,31,76,5,34;...
    0,21,13,48,11,59;...
    0,10,26,17,42,57;...
    0,7,45,-4,17,4;...
    0,8,5,-18,7,51];

obsv3 = [23,52,32.4619,-41,30,1.3241;...
    23,53,4.3124,-36,52,50.2637;...
    23,55,33.775,-31,4,37.9927;...
    0,0,50.7998,-23,21,5.6438;...
    0,10,31.5456,-12,18,30.1250;...
    0,27,58.4948,4,22,35.8242;...
    1,1,12.2853,28,16,30.6060;...
    2,11,7.8713,53,19,19.9219];

obsv_s = obsv1;

obsv = [obsv_s(:,1)*360/(24)+obsv_s(:,2)*360/(24*60)+obsv_s(:,3)*360/(24*3600),...
    obsv_s(:,4)+sign(obsv_s(:,4)).*(obsv_s(:,5)/60+obsv_s(:,6)/3600)];

t1 = [23,36,1.5;...
    23,37,0.8;...
    23,38,0;...
    23,39,1;...
    23,40,1.5;...
    23,40,58.5];

t2 = [23,16,10;...
    23,17,10.5;...
    23,18,11;...
    23,19,11;...
    23,20,12;...
    23,21,12];

t3 = [22,30,32;...
    22,31,32;...
    22,32,32;...
    22,33,32;...
    22,34,32;...
    22,35,32;...
    22,36,32;...
    22,37,32];

t = t1;

t_s = t(:,1)*60*60+t(:,2)*60+t(:,3);

obsv_select = [obsv(1,:);obsv(3,:);obsv(6,:)];
T_select = [t_s(1)-t_s(3);0;t_s(6)-t_s(3)];

lat = 39.00198;
lon = -76.9560;
R = zeros(3,size(t,1));
for i = [1:size(t,1)]
    ST_i = sidereal(2017,11,13,t(i,1),t(i,2),t(i,3),lon);
    ST_i = mod(ST_i,360)+0.5;
    Ri = SiteVec(lat,ST_i,0.055);
    R(:,i) = Ri;
end
R_select = [R(:,1),R(:,3),R(:,6)];

% lat = 40; alt = 1;
% obsv_select = [43.537,-8.7833;54.420,-12.074;64.318,-15.105];
% ST_1 = 44.506; ST_2 = 45; ST_3 = 45.499;
% R_select = [SiteVec(lat,ST_1,alt),SiteVec(lat,ST_2,alt),SiteVec(lat,ST_3,alt)];
% T_select = [-118.1,0,(237.58-118.10)];

r_set = AnglesGauss(obsv_select,R_select,T_select,mu);
r1 = r_set(:,1);
v1 = Gibbs(r_set,mu);
h1 = cross(r1,v1);
n = cross([1;0;0],h1);
o_e = orbital_elements(r1,v1,mu)

tol = 1E-13;
r_obsv = zeros(3,size(t,1));
r_prop = zeros(3,size(t,1));
v_prop = zeros(3,size(t,1));
u = zeros(size(t,1),1);
Period = 2*pi*sqrt((o_e(1)^3)/mu);
for i = [1:size(t,1)]
    t2 = t_s(i)-t_s(1);
%     if t2 < 0
%         t2 = double(floor(vpa(t2 + Period)));
%     end
    t_range = [1 t2];
    initial_conds = [r1',v1'];
    [t2,RV2] = Position_2BP(initial_conds,t_range,tol,mu);
    ri = RV2(size(t2,1),1:3)';
    r_prop(:,i) = ri;
    vi = RV2(size(t2,1),4:6)';
    v_prop(:,i) = vi;
    range_vec_i = ri-R(:,i); range_i = norm(range_vec_i);
    p_i = range_i*[cos(obsv(i,1)*pi/180)*cos(obsv(i,2)*pi/180);...
        cos(obsv(i,1)*pi/180)*sin(obsv(i,2)*pi/180);...
        sin(obsv(i,2)*pi/180)];
    r_i = p_i + R(:,i);
    r_obsv(:,i) = r_i;
    ui = acos(dot(n,r_i)/(norm(n)*norm(r_i)));
    u(i,1) = ui;
end

du = zeros(size(t,1)-1,1);
a = zeros(size(t,1)-1,1);
ddu = zeros(size(t,1)-1,1);
for i = [1:size(t,1)-1]
    du(i,1) = u(i+1)-u(i);
    dti = t_s(i+1)-t_s(i);
    a(i,1) = (((((2*pi/du(i,1))*dti)/(2*pi))^2)*mu)^(1/3);
    ddu(i,1) = -3*dti*sqrt(mu/(a(i,1)^3))/(2*a(i,1));
end

iterations = [1:10];
an = a(1,1)
for i = iterations
    uin = zeros(size(t,1)-1,1);
    A = zeros(size(t,1)-1,1);
    for j = [1:size(t,1)-1]
        dtj = t_s(j+1)-t_s(j);
        uin(j,1) = dtj*sqrt(mu/(an^3));
        A(j,1) = (-3*dtj/(2*an))*sqrt(mu/(an^3));
    end
    b = du-uin
    norm(b)
    dan = inv(A'*A)*A'*b;
    an = an+dan
end