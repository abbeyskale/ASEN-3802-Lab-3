%% ASEN 3111 - Computational Assigment 01 - Main
% Problem Statement

% Author: Maria Arrazola
% Collaborators: Abbey Skale, Anton Pylypchenko
% Date: March 27th 2026

clc; clear; close all;

%% Parameters 
% NACA 0021 Parameters 
m1= 0;
p1= 0;
t1= 21/100;

% NACA 2421 Parameters
m2= 2/100;
p2= 4/10;
t2= 21/100;

%% Call functions 
[x1,y1,~,~]=NACA_Airfoils(m1,p1,t1,c,N);
[x2,y2,xc,yc]=NACA_Airfoils(m2,p2,t2,c,N);

%% Plot 