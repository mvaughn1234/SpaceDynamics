clear all; close all; clc;

u = 3.986*(10^5);
tolerance = 10^-13;
trange = [1:30000];
input = [5203.12878457022 2539.18526782417 4387.98384076804 -5.73055171828814 1.23647597198147 6.07959326945700];

[t,RV] = Position_2BP(input,trange,tolerance,u);

for i = [1:30000]
    r(i) = norm(RV(i,1:3));
    v(i) = norm(RV(i,4:6));
    a(i) = u/(norm(RV(i,1:3))^2);
    
    sp(i) = (norm(RV(i,4:6))^2)/2 - u/norm(RV(i,1:3));
    h(i,1:3) = cross(RV(i,1:3),RV(i,4:6));
    hmag(i) = norm(cross(RV(i,1:3),RV(i,4:6)));
end

pos = RV(30000,1:3)
vel = RV(30000,4:6)

plot(t,sp);
figure;
subplot(2,2,1);
plot(t,hmag);
subplot(2,2,2);
plot(t,h(:,1));
subplot(2,2,3);
plot(t,h(:,2));
subplot(2,2,4);
plot(t,h(:,3));

