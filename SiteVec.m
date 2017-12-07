% Lat in degrees, sidereal time in degrees, altitude in km
function R = SiteVec(lat,sidereal_time,alt)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Re = 6378.1363;
ee = 0.081819;

xp = (Re/sqrt(1-(ee^2)*sin(lat*pi/180)^2)+alt)*cos(lat*pi/180);
zp = (Re*(1-ee^2)/sqrt(1-(ee^2)*sin(lat*pi/180)^2)+alt)*sin(lat*pi/180);

R = [xp*cos(sidereal_time*pi/180);xp*sin(sidereal_time*pi/180);zp];
end

