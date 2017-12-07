clear all; close all; clc;

rA = [9222.33068385812 1408.28350774150 813.072862289868];
vA = [-1.29171353392020 6.34421760239082 3.66283574053790];
r = norm(rA), v = norm(vA), u = 3.986*10^5;
vesc = sqrt(2*u/r)

figure; hold on;
plot_v(r,v);
syms g(x);
hold off;

function nothing = plot_v(r,v)
u = 3.986*10^5;
syms f(x);
E = (v^2)/2 - u/r
a = -u/(2*E)
e = 1 - r/a
b = sqrt((a^2)*(1-e^2))

f(x) = sqrt(b^2-(b*(x-a)/a)^2);
fplot(-f(x),[0 2*a]);
fplot(f(x),[0 2*a]);
end