% First convert from elements to PQW then from PQW to ECI then from ECI To
% ECEF
a = 8350;
e = 0.1;
i = 41;
w = 295;
omega = 302;
nu = 3.986*10^5;
M = 22.2;
elem = [a,e,i,omega,w,M];
eciRV = elementsToECI(elem,nu)
d = 5007+10/(60*24)+22/24;
th = 280.4606+360.9856473*d
thrad = th*pi/180;
R3th = [[1 0 0];[0 cos(thrad) sin(thrad)];[0 -sin(thrad) cos(thrad)]];
ecefR = R3th*eciRV(1:3)

