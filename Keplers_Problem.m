clear all; close all; clc; format compact; format long;
case1 = [15000,0.3,60,30,0,120];
case2 = [21000,0.7,20,0,0,35];
Uearth = 3.986*10^5;
eci1 = elementsToECI(case1,Uearth)
eci2 = elementsToECI(case2,Uearth)
t = 1*60*60;
[r2_1,v2_1,t1] = calc_keplers(case1,eci1(1:3),eci1(4:6),xiGuess(case1,t,Uearth),Uearth,t,1*10^-10)
[r2_2,v2_2,t2] = calc_keplers(case2,eci2(1:3),eci2(4:6),xiGuess(case2,t,Uearth),Uearth,t,1*10^-10)
[tcheck1,RVcheck1] = Position_2BP(eci1,[0,t],1*10^-13,Uearth);