clear all; close all; clc;
u_e = 398600.441800;

A30 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\A30 Inertial Position Velocity.csv');
A60 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\A60 Inertial Position Velocity.csv');
A90 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\A90 Inertial Position Velocity.csv');
A120 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\A120 Inertial Position Velocity.csv');
A150 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\A150 Inertial Position Velocity.csv');
B30 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\B30 Inertial Position Velocity.csv');
B60 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\B60 Inertial Position Velocity.csv');
B90 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\B90 Inertial Position Velocity.csv');
B120 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\B120 Inertial Position Velocity.csv');
B150 = dlmread('C:\Users\mvaug\Documents\STK 11 (x64)\STK_Project\B150 Inertial Position Velocity.csv');
all_plots = [A30 A60 A90 A120 A150 B30 B60 B90 B120 B150];
for plot_num = [1:1:10]
    left_column = (plot_num-1)*7+1;
    right_column = plot_num*7;
    current_plot = all_plots(:,left_column:right_column);
    for i = [1:1:5761]
       t(i) = current_plot(i,1);
       eles(i,1:6) = orbital_elements(current_plot(i,2:4)',current_plot(i,5:7)',u_e);
       if(i > 1)
           prev = i-1;
           regression(prev,1) = eles(i,3);
           regression(prev,2) = ((eles(i,4)-eles(prev,4))/60)*(24*60*60);
           rotation(prev,1) = eles(i,3);
           rotation(prev,2) = ((eles(i,5)-eles(prev,5))/60)*(24*60*60);
       end
    end
%     if(plot_num == 1)
        figure; hold on;
        subplot(2,3,1);
        plot(t,eles(:,1))
        title('a');
        xlabel('Time [s]');
        ylabel('[km]');

        subplot(2,3,2);
        plot(t,eles(:,2))
        title('e');
        xlabel('Time [s]');

        subplot(2,3,3);
        plot(t,eles(:,3))
        title('i');
        xlabel('Time [s]');
        ylabel('[degrees]');

        subplot(2,3,4);
        plot(t,eles(:,4))
        title('Omega');
        xlabel('Time [s]');
        ylabel('[degrees]');

        subplot(2,3,5);
        plot(t,eles(:,5))
        title('w');
        xlabel('Time [s]');
        ylabel('[degrees]');

        subplot(2,3,6);
        plot(t,eles(:,6))
        title('Nu');
        xlabel('Time [s]');
        ylabel('[degrees]');
        hold off;
%     end
    if(plot_num==1)
        avg_regression = mean(regression(:,2));
        avg_reg_total = [30,avg_regression];
    elseif(plot_num==2)
        avg_regression = mean(regression(:,2));
        avg_reg_total = [avg_reg_total;60,avg_regression];
    elseif(plot_num==3)
        avg_regression = mean(regression(:,2));
        avg_reg_total = [avg_reg_total;90,avg_regression];
    elseif(plot_num==4)
        avg_regression = mean(regression(:,2));
        avg_reg_total = [avg_reg_total;120,avg_regression];
    elseif(plot_num==5)
        avg_regression = mean(regression(:,2));
        avg_reg_total = [avg_reg_total;150,avg_regression];
    elseif(plot_num==6)
        avg_regression = mean(regression(:,2));
        avg_rotation = mean(rotation(:,2));
        avg_reg_total2 = [30,avg_regression];
        avg_rot_total2 = [30,avg_rotation];
    elseif(plot_num==7)
        avg_regression = mean(regression(:,2));
        avg_rotation = mean(rotation(:,2));
        avg_reg_total2 = [avg_reg_total2;60,avg_regression];
        avg_rot_total2 = [avg_rot_total2;60,avg_rotation];
    elseif(plot_num==8)
        avg_regression = mean(regression(:,2));
        avg_rotation = mean(rotation(:,2));
        avg_reg_total2 = [avg_reg_total2;90,avg_regression];
        avg_rot_total2 = [avg_rot_total2;90,avg_rotation];
    elseif(plot_num==9)
        avg_regression = mean(regression(:,2));
        avg_rotation = mean(rotation(:,2));
        avg_reg_total2 = [avg_reg_total2;120,avg_regression];
        avg_rot_total2 = [avg_rot_total2;120,avg_rotation];
    elseif(plot_num==10)
        avg_regression = mean(regression(:,2));
        avg_rotation = mean(rotation(:,2));
        avg_reg_total2 = [avg_reg_total2;150,avg_regression];
        avg_rot_total2 = [avg_rot_total2;150,avg_rotation];
    end
end
figure; hold on;
plot(avg_reg_total(:,1),avg_reg_total(:,2));
plot(avg_reg_total2(:,1),avg_reg_total2(:,2));
hold off;
xlabel('Inclination [deg]'); ylabel('Nodal Regression [deg/day]');
legend('A','B','Location','NorthWest');
figure;
plot(avg_rot_total2(:,1),avg_rot_total2(:,2));
xlabel('Inclination [deg]'); ylabel('Apsidal Rotation [deg/day]');
legend('B','Location','NorthWest');