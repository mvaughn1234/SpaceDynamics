%% ENAE 441 Final Project
% Mathew Vaughn
clear all; close all; clc;

Re = 6371; % km
mu = 3.986E5;

obsv = [21*360/24+14*360/(24*60),-28;...
    21*360/24+15*360/(24*60),-15;...
    21*360/24+15*360/(24*60),8+30*360/(24*60*60);...
    21*360/24+15*360/(24*60),47+15*360/(24*60*60);...
    21*360/24+10*360/(24*60),81;...
    9*360/24+30/(24*60),81];

t = [23,36,1.5;...
    23,37,0.8;...
    23,38,0;...
    23,39,1;...
    23,40,1.5;...
    23,40,58.5];

t_s = t(:,1)*24*60+t(:,2)*60+t(:,3);

obsv_select = [obsv(1,:);obsv(3,:);obsv(6,:)];
T_select = [t_s(1)-t_s(3);0;t_s(6)-t_s(3)];

lat = 39.00198;
lon = -76.9560;
R = zeros(3,size(t,1));
for i = [1:size(t,1)]
    ST_i = sidereal(2017,11,13,t(i,1),t(i,2),t(i,3),lon);
    Ri = SiteVec(lat,ST_i,0.055);
    R(:,i) = Ri;
end
R_select = [R(:,1),R(:,3),R(:,6)];

r_set = AnglesGauss(obsv_select,R_select,T_select,mu);
v2 = Gibbs(r_set,mu);
o_e = orbital_elements(r2,v2,mu);

tol = 1E-13;
r_obsv = zeros(3,size(t,1));
for i = [1:size(t,1)]
    t_range = [t_s(3):t_s(i)];
    [t2,RV2] = Position_2BP([r2',v2'],t_range,tol,mu)
    ri = RV2(size(t2,1),1:3);
    vi = RV2(size(t2,1),4:6);
    range_vec_i = ri-R(:,i); range_i = norm(range_vec_i);
    pi = range*[cos(obsv(i,1)*pi/180)*cos(obsv(i,2)*pi/180);...
        cos(obsv(i,1)*pi/180)*sin(obsv(i,2)*pi/180);...
        sin(obsv(i,2)*pi/180)];
    r_i = pi + R(:,i);
    r_obsv(:,i) = r_i;
    
    
    
    
    