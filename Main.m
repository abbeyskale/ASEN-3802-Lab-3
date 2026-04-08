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


%% -------- Task 3 ------------

% For symmetric airfoils, alpha_L=0 = 0 deg and:
% c_l = 2*pi*alpha, where alpha must be in radians
alphaTheory_deg = -15:0.1:15;
clTheory = 2*pi*deg2rad(alphaTheory_deg);
slopeTheory_deg = 2*pi*(pi/180);

%% Experimental data: NACA 0006
exp0006 = [
  -10.2824   -0.871349
   -8.00549  -0.794991
   -6.33373  -0.665445
   -4.96550  -0.543543
   -3.90121  -0.444504
   -2.37988  -0.269323
   -0.706679 -0.0865151
    0.814653  0.0886658
    2.94345   0.294354
    4.16043   0.431455
    5.83301   0.591436
    7.05020   0.736146
    8.41802   0.842831
];

alpha0006_deg = exp0006(:,1);
cl0006_exp = exp0006(:,2);

% Selected linear region for slope calculation
linearMask0006 = (alpha0006_deg >= -6) & (alpha0006_deg <= 6);
alpha0006_lin = alpha0006_deg(linearMask0006);
cl0006_lin = cl0006_exp(linearMask0006);

% Linear fit: c_l = m*alpha + b
fit0006 = polyfit(alpha0006_lin, cl0006_lin, 1);
slope0006_exp = fit0006(1);
alphaL0_0006_exp = -fit0006(2)/fit0006(1);

%% Experimental data: NACA 0012
exp0012 = [
  -15.7811   -1.55776
  -14.2488   -1.47278
  -13.1343   -1.37389
  -11.0448   -1.17616
  -10.0697   -1.08435
   -8.11940  -0.851509
   -6.16915  -0.667899
   -5.05473  -0.526806
   -2.82587  -0.307924
   -0.875622 -0.110247
    1.07463   0.115564
    3.02488   0.341376
    5.25373   0.532124
    7.90050   0.877769
   10.1294    1.08258
   12.0796    1.24509
   14.0299    1.41464
   15.5622    1.54885
];

alpha0012_deg = exp0012(:,1);
cl0012_exp = exp0012(:,2);

% Selected linear region for slope calculation
linearMask0012 = (alpha0012_deg >= -8) & (alpha0012_deg <= 5);
alpha0012_lin = alpha0012_deg(linearMask0012);
cl0012_lin = cl0012_exp(linearMask0012);

% Linear fit: c_l = m*alpha + b
fit0012 = polyfit(alpha0012_lin, cl0012_lin, 1);
slope0012_exp = fit0012(1);
alphaL0_0012_exp = -fit0012(2)/fit0012(1);

%% Print experimental and thin airfoil theory results to command window
fprintf('Experimental lift slope for NACA 0006 = %.4f per degree\n', slope0006_exp);
fprintf('Experimental zero-lift angle for NACA 0006 = %.4f deg\n', alphaL0_0006_exp);
fprintf('\n');

fprintf('Experimental lift slope for NACA 0012 = %.4f per degree\n', slope0012_exp);
fprintf('Experimental zero-lift angle for NACA 0012 = %.4f deg\n', alphaL0_0012_exp);
fprintf('\n');

fprintf('Thin airfoil theory lift slope = %.4f per degree\n', slopeTheory_deg);
fprintf('Thin airfoil theory zero-lift angle for symmetric airfoils = 0.0000 deg\n');
fprintf('\n');

%% Vortex panel method setup

% Use the total number of panels determined in Task 2
c = 1;
N_panels = 40;
alphaVP_deg = -15:1:15;

% Symmetric airfoil parameters
m = 0;
p = 0;

t0006 = 0.06;
t0012 = 0.12;
t0018 = 0.18;

%% Generate airfoil boundary points for vortex panel method
[xb0006, yb0006, ~, ~] = NACA_Airfoils(m, p, t0006, c, N_panels);
[xb0012, yb0012, ~, ~] = NACA_Airfoils(m, p, t0012, c, N_panels);
[xb0018, yb0018, ~, ~] = NACA_Airfoils(m, p, t0018, c, N_panels);

%% Run vortex panel method over angle of attack range
cl0006_vp = zeros(size(alphaVP_deg));
cl0012_vp = zeros(size(alphaVP_deg));
cl0018_vp = zeros(size(alphaVP_deg));

for i = 1:length(alphaVP_deg)
    alpha_i = alphaVP_deg(i);

    cl0006_vp(i) = Vortex_Panel_2(xb0006, yb0006, alpha_i);
    cl0012_vp(i) = Vortex_Panel_2(xb0012, yb0012, alpha_i);
    cl0018_vp(i) = Vortex_Panel_2(xb0018, yb0018, alpha_i);
end

%% Compute vortex panel lift slopes and zero-lift angles
linearMaskVP = (alphaVP_deg >= -6) & (alphaVP_deg <= 6);

fit0006_vp = polyfit(alphaVP_deg(linearMaskVP), cl0006_vp(linearMaskVP), 1);
fit0012_vp = polyfit(alphaVP_deg(linearMaskVP), cl0012_vp(linearMaskVP), 1);
fit0018_vp = polyfit(alphaVP_deg(linearMaskVP), cl0018_vp(linearMaskVP), 1);

slope0006_vp = fit0006_vp(1);
slope0012_vp = fit0012_vp(1);
slope0018_vp = fit0018_vp(1);

alphaL0_0006_vp = -fit0006_vp(2)/fit0006_vp(1);
alphaL0_0012_vp = -fit0012_vp(2)/fit0012_vp(1);
alphaL0_0018_vp = -fit0018_vp(2)/fit0018_vp(1);

%% Print vortex panel results to command window
fprintf('Vortex panel lift slope for NACA 0006 = %.4f per degree\n', slope0006_vp);
fprintf('Vortex panel zero-lift angle for NACA 0006 = %.4f deg\n', alphaL0_0006_vp);
fprintf('\n');

fprintf('Vortex panel lift slope for NACA 0012 = %.4f per degree\n', slope0012_vp);
fprintf('Vortex panel zero-lift angle for NACA 0012 = %.4f deg\n', alphaL0_0012_vp);
fprintf('\n');

fprintf('Vortex panel lift slope for NACA 0018 = %.4f per degree\n', slope0018_vp);
fprintf('Vortex panel zero-lift angle for NACA 0018 = %.4f deg\n', alphaL0_0018_vp);

%% Plot 

%Task 1 Plots 
figure(1);
legend('NACA 2421', 'Camber Line');
plot(x1,y1, 'LineStyle','-');
hold on;
scatter(x1,y1,'filled');
title('NACA 0021');
yline(0,'Color','r','LineStyle','--'); % camber line
axis equal; 
legend('Airfoil Shape','Panel Separations','Camber Line');
grid on;
print('Task1NACA0021','-dpng','-r300');
hold off

figure(2);
plot(x2,y2,'LineStyle','-'); 
hold on;
scatter(x2,y2,'filled');
plot(x_cam,yc,'Color','r','LineStyle','--'); % camber line
title('NACA 2421 with Camber Line');
axis equal; 
legend('Airfoil Shape','Panel Separations','Camber Line');
grid on;
print('Task1NACA2421','-dpng','-r300')
hold off

%Task 2 Plots

figure(3);
plot(N_val,CL_val,'-o');
hold on;
yline(CL3,'Linestyle','--');
xline(N_min,'Linestyle','--');
xlabel('Number of Panels');
ylabel('Cl');
title('Convergence of cl with Number of Panels');
legend('Computed Cl', 'Exact Cl','1% error', 'Location','Best');
grid on;
hold off

% FOR NACA 0012 Geometry 
figure(4);
plot(x3,y3,'LineStyle','-');
hold on;
scatter(x3,y3,'filled');
plot(x_cam,yc,'Color','r','LineStyle','--'); 
title('NACA 0012');
axis equal;
legend('Airfoil Shape','Panel Separations','Camber Line');
grid on;
hold off

%% Task 3 Plots
figure(4)
plot(alpha0006_deg, cl0006_exp, 'r-', 'LineWidth', 1.5, 'MarkerSize', 5)
hold on
plot(alpha0012_deg, cl0012_exp, 'b-', 'LineWidth', 1.5, 'MarkerSize', 5)
plot(alphaTheory_deg, clTheory, 'k-', 'LineWidth', 1.8)

plot(alphaVP_deg, cl0006_vp, 'r--', 'LineWidth', 1.8)
plot(alphaVP_deg, cl0012_vp, 'b--', 'LineWidth', 1.8)
plot(alphaVP_deg, cl0018_vp, 'g--', 'LineWidth', 1.8)

xlabel('\alpha (deg)')
ylabel('c_l')
title('Experimental, Thin Airfoil Theory, and Vortex Panel Lift Curves')
legend('NACA 0006 Experimental', ...
       'NACA 0012 Experimental', ...
       'Thin Airfoil Theory', ...
       'NACA 0006 Vortex Panel', ...
       'NACA 0012 Vortex Panel', ...
       'NACA 0018 Vortex Panel', ...
       'Location', 'best')
grid on
xlim([-16 16])
hold off

%% NACA Airfoils function
function [x_b, y_b, x, yc] = NACA_Airfoils(m,p,t,c,N)

%cosine spacing 
theta=linspace(0,pi,N);
x=(1/2)*(1-cos(theta)); 

%Thicknesss distribution 
yt= (c*t/0.2)*((0.2969*sqrt(x/c))- (0.126*(x/c))-(0.3516*(x/c).^2)+(0.2843*(x/c).^3)-(0.1036*(x/c).^4));

% camber line and slope as vectors
if p==0
yc= zeros(size(x));
dycdx= zeros(size(x));
else
   if x < p*c
    yc = (m.*x./p.^2)*(2.*p-x./c);
    dycdx = 2.*m.*x./p - 2*m./(p.^2.*c);
   else
    yc = m.*(c-x)./(1-p).^2 .* (1+ x./c - 2*p);
    dycdx= 2*m.*x/((1-p)^2*c) - 2*m*p/(1-p)^2;
   end 
end

% Angle
Xi= atan(dycdx);

% Upper Surface
xU= x-yt .*sin(Xi);
yU= yc+yt .*cos(Xi);

%Lower Surface
xL= x+yt .*sin(Xi);
yL= yc-yt.*cos(Xi);

%% Combine surfaces into one (clockwise)
% Upper surface tailing edge to leading edge
xU= flip(xU);
yU= flip(yU);

% Lower surface leading edge to trailing edge
xL= xL(2:end);   %skip 1st point of the lower surface so doesn't duplicate Leaing edge
yL= yL(2:end);




% Final boundary points
x_b= [xU,xL];
y_b= [yU,yL];
y_b = flip(y_b);
end

