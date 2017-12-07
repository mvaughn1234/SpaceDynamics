function[r,v]=findEarth(JD)
%%%%gives the inertial position and velocity of Earth at a given Julian day
%%%% created by Dr. Christine Hartzell (2010)
%%%% Distributed to ENAE 404 for use in project assignment
%%%% update with your own oe2cart code

TDB=(JD-2451545)/36525;
AU=1.4959787E8; %km
mu=132712440018; %km^3/s^2 (this is for the sun)
a=1.000001018*AU;
e=0.01670862- 0.000042037*TDB- 0.0000001236*TDB^2+0.00000000004*TDB^3;
i=0.0130546*TDB- 0.00000931*TDB^2 - 0.000000034*TDB^3;%deg
Omega=174.873174-0.2410908*TDB+0.00004067*TDB^2-0.000001327*TDB^3;%deg
wtilde=102.937348+0.3225557*TDB+0.00015026*TDB^2+0.000000478*TDB^3; %deg
L=100.466449+35999.3728519*TDB- 0.00000568*TDB^2; %deg

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