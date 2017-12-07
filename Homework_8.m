clear all; close all; clc;
I = [10,0,0;0,20,0;0,0,30];
w = [5;0;25];
t = [0 100];
tolerance = 1*10^-13;
[t,wnew] = Prop_Ang_Mom(I,w,t,tolerance);
figure; hold on;
plot(t,wnew(:,1),'black')
plot(t,wnew(:,2),'black--')
plot(t,wnew(:,3),'black-.')
hold off;
length = size(t,1)
for i = [1:1:length]
    T(i) = 0.5*I(1,1)*wnew(i,1)^2+0.5*I(2,2)*wnew(i,2)^2+0.5*I(3,3)*wnew(i,3)^2;
    H(i) = sqrt(I(1,1)^2*wnew(i,1)^2+I(2,2)^2*wnew(i,2)^2+I(3,3)^2*wnew(i,3)^2);
end
figure; hold on;
subplot(2,1,1)
plot(t,T,'black')
subplot(2,1,2)
plot(t,H,'black')
hold off;
