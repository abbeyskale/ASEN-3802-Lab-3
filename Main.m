%% ASEN 3802 - Computational Assigment 01 - Main
% This code uses the vortex panel method to plot the x and y coordinates of
% symmetrical and cambered airfoils

% Author: Maria Arrazola
% Collaborators: Abbey Skale, Anton Pylypchenko
% Date: March 27th 2026

clc; clear; close all;

%% Parameters 
c=1;
N=50;
t=21/100;

% Flow conditions
alpha = 12; % degrees

%Task 1 Parameters 
% NACA 0021 Parameters 
m1= 0;
p1= 0;

% NACA 2421 Parameters
m2= 2/100;
p2= 4/10;

% Task 2 NACA 0012 Parameters
m3=0;
p3=0;
N3=500;
t3= 12/100;

%% Call functions 
%Generating airfoils
[x1,y1,~,~]=NACA_Airfoils(m1,p1,t,c,N);
[x2,y2,x_cam,yc]=NACA_Airfoils(m2,p2,t,c,N);
[x3,y3,~,~]=NACA_Airfoils(m3,p3,t3,c,N3);

%Running vortex panel method 
CL1 = Vortex_Panel_2(x1,y1,alpha);
CL2 = Vortex_Panel_2(x2,y2,alpha);
CL3 = Vortex_Panel_2(x3,y3,alpha);

%% Displaying CL for Airfoils
%fprintf('CL for NACA 0021: %.4f\n', CL1);
%fprintf('CL for NACA 2421: %.4f\n', CL2);
fprintf('CL for NACA 0012: %.4f\n using N= %.0f\n' , CL3,N3);


%% Testing N values for convergence
N_val=[4,8,12,16,20,24,28,32,36,40,44,48,52];
CL_val= zeros(size(N_val));
error=zeros(size(N_val));

for i=1:length(N_val)
    N_tot= N_val(i);
    [xb,yb,~,~]= NACA_Airfoils(m3,p3,t3,c,N_tot);
    CL_val(i)= Vortex_Panel_2(xb,yb,alpha);
    error(i)= abs(CL_val(i)- CL3)/CL3*100;
end

% Find minimum N for less than 1% relative error
index= find(error < 1,1,'first');
N_min= N_val(index);

fprintf('Minimum N for less than 1 percent relative error is %d panels\n',N_min);
fprintf('CL at N= %d is %.4f, relative error is %.2f%%\n', N_min,CL_val(index),error(index));


%% Plot 

%Task 1 Plots 
figure;
legend('NACA 2421', 'Camber Line');
plot(x1,y1, 'LineStyle','-');
hold on;
scatter(x1,y1,'filled');
title('NACA 0021');
yline(0,'Color','r','LineStyle','--'); % camber line
axis equal; 
legend('Airfoil Shape','Panel Separations','Camber Line');
grid on;
print('Task1NACA0021','-dpng','-r300')

figure;
plot(x2,y2,'LineStyle','-'); 
hold on;
scatter(x2,y2,'filled');
plot(x_cam,yc,'Color','r','LineStyle','--'); % camber line
title('NACA 2421 with Camber Line');
axis equal; 
legend('Airfoil Shape','Panel Separations','Camber Line');
grid on;
print('Task1NACA2421','-dpng','-r300')

%Task 2 Plots
figure;
plot(x3,y3,'LineStyle','-');
hold on;
scatter(x3,y3,'filled');
plot(x_cam,yc,'Color','r','LineStyle','--'); 
title('NACA 0012');
axis equal;
legend('Airfoil Shape','Panel Separations','Camber Line');
grid on;




