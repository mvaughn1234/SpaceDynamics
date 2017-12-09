function [ ThL ] = sidereal(year, month, day, hour, minute, second, Lon)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
jd = ymdhms2jd(year,month,day,hour,minute,second);
djd = jd-ymdhms2jd(2000,1,1,12,0,0);
jct = djd/36525;

ThL = 100.4606184+36000.7700*jct+0.000387933*jct^2-(2.583E-8)*jct^3 + ...
    360.985647*((hour+minute/60+second/3600)/24) + Lon;

end

