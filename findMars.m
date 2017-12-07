function[r,v]=findMars(JD)
%%%%gives the inertial position and velocity of Earth at a given Julian day
%%%% created by Dr. Christine Hartzell (2010)
%%%% Distributed to ENAE 404 for use in project assignment
%%%% update with your own oe2cart code
TDB=(JD-2451545)/36525;
AU=1.4959787E8; %km
mu=132712440018; %km^3/s^2 (this is for the sun)

a=1.523679342*AU;
e=0.09340062+0.000090483*TDB- 0.0000000806*TDB^2 - 0.00000000035*TDB^3;
i= 	 1.849726 - 0.0081479*TDB - 0.00002255*TDB^2- 0.000000027*TDB^3;%deg
Omega=49.558093- 0.2949846*TDB- 0.00063993*TDB^2- 0.000002143*TDB^3;%deg
wtilde=336.060234+0.4438898*TDB - 0.00017321*TDB^2+0.000000300*TDB^3; %deg
L=355.433275+19140.2993313*TDB+0.00000261*TDB^2- 0.000000003*TDB^3; %deg

M=L-wtilde;
M=rem(M,360);%deg
w=wtilde-Omega;

E=findE(M*pi/180,e);%radians
E_check=rem(E,2*pi);
nu=acos((cos(E)-e)/(1-e*cos(E)));%rad
%%%Check E and nu should be in the same half plane
if (E_check>pi &&nu<pi) |(E_check<pi && nu>pi)
    nu=2*pi-nu;
end

[r,v]=oe2cart(a,e,i,Omega,w,nu*180/pi);