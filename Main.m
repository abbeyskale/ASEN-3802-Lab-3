%% ASEN 3802 - Computational Assigment 01 - Main
% This code uses the vortex panel method to plot the x and y coordinates of
% symmetrical and cambered airfoils

% Author: Maria Arrazola
% Collaborators: Abbey Skale
% Date: March 27th 2026

clc; clear; close all;

%% Parameters 
c=1;
N=50;
t=21/100;

% Flow conditions
alpha = 12; % degrees

% NACA 0021 Parameters 
m1= 0;
p1= 0;

% NACA 2421 Parameters
m2= 2/100;
p2= 4/10;

%% Call functions 
%Generating airfoils
[x1,y1,~,~]=NACA_Airfoils(m1,p1,t,c,N);
[x2,y2,x_cam,yc]=NACA_Airfoils(m2,p2,t,c,N);

%Running vortex panel method 
CL1 = Vortex_Panel_2(x1,y1,alpha);
CL2 = Vortex_Panel_2(x2,y2,alpha);

%% Displaying CL for Airfoils 
fprintf('CL for NACA 0021: %.4f\n', CL1);
fprintf('CL for NACA 2421: %.4f\n', CL2);

%% Plot 
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

figure;
plot(x2,y2,'LineStyle','-'); 
hold on;
scatter(x2,y2,'filled');
plot(x_cam,yc,'Color','r','LineStyle','--'); % camber line
title('NACA 2421 with Camber Line');
axis equal; 
legend('Airfoil Shape','Panel Separations','Camber Line');
grid on;



