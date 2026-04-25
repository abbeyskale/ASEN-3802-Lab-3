%% ASEN 3802 - Computational Assigment 01 - Main
% This code uses the vortex panel method to plot the x and y coordinates of
% symmetrical and cambered airfoils

% Author: Maria Arrazola
% Collaborators: Abbey Skale, Anton Pylypchenko
% Date: April 15th 2026

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
fprintf('CL at N= %d is %.4f, relative error is %.2f%%\n\n', N_min,CL_val(index),error(index));


%% -------- Task 3 ------------

% For symmetric airfoils, alpha_L=0 = 0 deg and 
% c_l = 2*pi*alpha, where alpha is in radians
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
N_panels = 36;
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

%% ---- Task 4 -----

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


exp2412 = [
 -12.3333     -1.00098
 -10.1667    -0.808754
 -8.33333    -0.625012
 -6.33333    -0.416259
 -4.16667    -0.190801
       -2    0.0429649
-0.166667     0.243323
  2.16667     0.427324
  4.16667      0.62777
  6.16667     0.836523
        8      1.02026
  10.3333      1.16272];

alpha2412_deg = exp2412(:,1);
cl2412_exp = exp2412(:,2);

% Selected linear region for slope calculation
linearMask2412 = (alpha2412_deg >= -8) & (alpha2412_deg <= 6);
alpha2412_lin = alpha2412_deg(linearMask2412);
cl2412_lin = cl2412_exp(linearMask2412);

% Linear fit: c_l = m*alpha + b
fit2412 = polyfit(alpha2412_lin, cl2412_lin, 1);
slope2412_exp = fit2412(1);
alphaL0_2412_exp = -fit2412(2)/fit2412(1);

exp4412 = [
-12.0235    -0.761387
-10.0614    -0.670427
-7.76352    -0.464257
-5.79063    -0.249946
-3.81701   -0.0274128
-1.84412     0.186897
0.454483     0.401291
 2.10095     0.607294
 4.39956     0.821688
 6.69816      1.03608
 8.50461      1.20923
 10.6368      1.38247];

alpha4412_deg = exp4412(:,1);
cl4412_exp = exp4412(:,2);

% Selected linear region for slope calculation
linearMask4412 = (alpha4412_deg >= -8) & (alpha4412_deg <= 6);
alpha4412_lin = alpha4412_deg(linearMask4412);
cl4412_lin = cl4412_exp(linearMask4412);

% Linear fit: c_l = m*alpha + b
fit4412 = polyfit(alpha4412_lin, cl4412_lin, 1);
slope4412_exp = fit4412(1);
alphaL0_4412_exp = -fit4412(2)/fit4412(1);

%% Minimum Panel count
N_panels = 36;

% Angle ranges
alphaVP_deg = -15:1:15;
alphaTheory_deg = -15:0.1:15;

% Airfoil definitions
c = 1;
t = 0.12;

m0012 = 0.00;  p0012 = 0.0;
m2412 = 0.02;  p2412 = 0.4;
m4412 = 0.04;  p4412 = 0.4;

%% Thin airfoil theory

% Lift slope (same for all airfoils)
slopeTheory_deg = 2*pi*(pi/180);

% Zero-lift angle
alphaL0_0012_theory = 0;

alphaL0_2412_theory = -2*m2412/p2412*(1 - p2412);
alphaL0_4412_theory = -2*m4412/p4412*(1 - p4412);

% Convert to degrees
alphaL0_2412_theory = rad2deg(alphaL0_2412_theory);
alphaL0_4412_theory = rad2deg(alphaL0_4412_theory);

% Lift curves
cl0012_theory = 2*pi*deg2rad(alphaTheory_deg - alphaL0_0012_theory);
cl2412_theory = 2*pi*deg2rad(alphaTheory_deg - alphaL0_2412_theory);
cl4412_theory = 2*pi*deg2rad(alphaTheory_deg - alphaL0_4412_theory);

%% Vortex panel setup

[xb0012, yb0012, ~, ~] = NACA_Airfoils(m0012,p0012,t,c,N_panels);
[xb2412, yb2412, ~, ~] = NACA_Airfoils(m2412,p2412,t,c,N_panels);
[xb4412, yb4412, ~, ~] = NACA_Airfoils(m4412,p4412,t,c,N_panels);

cl0012_vp = zeros(size(alphaVP_deg));
cl2412_vp = zeros(size(alphaVP_deg));
cl4412_vp = zeros(size(alphaVP_deg));

for i = 1:length(alphaVP_deg)
    alpha_i = alphaVP_deg(i);
    cl0012_vp(i) = Vortex_Panel_2(xb0012,yb0012,alpha_i);
    cl2412_vp(i) = Vortex_Panel_2(xb2412,yb2412,alpha_i);
    cl4412_vp(i) = Vortex_Panel_2(xb4412,yb4412,alpha_i);
end

%% Vortex panel lift slope and zero-lift angle

linearMaskVP = (alphaVP_deg >= -6) & (alphaVP_deg <= 6);

fit0012_vp = polyfit(alphaVP_deg(linearMaskVP), cl0012_vp(linearMaskVP), 1);
fit2412_vp = polyfit(alphaVP_deg(linearMaskVP), cl2412_vp(linearMaskVP), 1);
fit4412_vp = polyfit(alphaVP_deg(linearMaskVP), cl4412_vp(linearMaskVP), 1);

slope0012_vp = fit0012_vp(1);
slope2412_vp = fit2412_vp(1);
slope4412_vp = fit4412_vp(1);

alphaL0_0012_vp = -fit0012_vp(2)/fit0012_vp(1);
alphaL0_2412_vp = -fit2412_vp(2)/fit2412_vp(1);
alphaL0_4412_vp = -fit4412_vp(2)/fit4412_vp(1);

%% Print results to command window

fprintf('\n');
fprintf('TASK 4 RESULTS\n');
fprintf('Using N = %d total panels from Task 2\n\n', N_panels);

fprintf('NACA 0012\n');
fprintf('Experimental lift slope = %.4f per degree\n', slope0012_exp);
fprintf('Experimental zero-lift angle = %.4f deg\n', alphaL0_0012_exp);
fprintf('Thin airfoil lift slope = %.4f per degree\n', slopeTheory_deg);
fprintf('Thin airfoil zero-lift angle = %.4f deg\n', alphaL0_0012_theory);
fprintf('Vortex panel lift slope = %.4f per degree\n', slope0012_vp); % should be 0.1189
fprintf('Vortex panel zero-lift angle = %.4f deg\n\n', alphaL0_0012_vp); 

fprintf('NACA 2412\n');
fprintf('Experimental lift slope = %.4f per degree\n', slope2412_exp);
fprintf('Experimental zero-lift angle = %.4f deg\n', alphaL0_2412_exp);
fprintf('Thin airfoil lift slope = %.4f per degree\n', slopeTheory_deg);
fprintf('Thin airfoil zero-lift angle = %.4f deg\n', alphaL0_2412_theory);
fprintf('Vortex panel lift slope = %.4f per degree\n', slope2412_vp);
fprintf('Vortex panel zero-lift angle = %.4f deg\n\n', alphaL0_2412_vp);

fprintf('NACA 4412\n');
fprintf('Experimental lift slope = %.4f per degree\n', slope4412_exp);
fprintf('Experimental zero-lift angle = %.4f deg\n', alphaL0_4412_exp);
fprintf('Thin airfoil lift slope = %.4f per degree\n', slopeTheory_deg);
fprintf('Thin airfoil zero-lift angle = %.4f deg\n', alphaL0_4412_theory);
fprintf('Vortex panel lift slope = %.4f per degree\n', slope4412_vp);
fprintf('Vortex panel zero-lift angle = %.4f deg\n\n', alphaL0_4412_vp);

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

% Task 4 Plot

figure(5);
plot(alpha0012_deg, cl0012_exp, 'ko', 'LineWidth', 1.5, 'MarkerSize', 6)
hold on
plot(alpha2412_deg, cl2412_exp, 'bo', 'LineWidth', 1.5, 'MarkerSize', 6)
plot(alpha4412_deg, cl4412_exp, 'ro', 'LineWidth', 1.5, 'MarkerSize', 6)

plot(alphaTheory_deg, cl0012_theory, 'k-', 'LineWidth', 1.5)
plot(alphaTheory_deg, cl2412_theory, 'b-', 'LineWidth', 1.5)
plot(alphaTheory_deg, cl4412_theory, 'r-', 'LineWidth', 1.5)

plot(alphaVP_deg, cl0012_vp, 'k--', 'LineWidth', 1.8)
plot(alphaVP_deg, cl2412_vp, 'b--', 'LineWidth', 1.8)
plot(alphaVP_deg, cl4412_vp, 'r--', 'LineWidth', 1.8)

xlabel('\alpha (deg)')
ylabel('c_l')
title('Effect of Airfoil Camber on Lift')
lgnd = legend('NACA 0012 Experimental', ...
       'NACA 2412 Experimental', ...
       'NACA 4412 Experimental', ...
       'NACA 0012 Thin Airfoil Theory', ...
       'NACA 2412 Thin Airfoil Theory', ...
       'NACA 4412 Thin Airfoil Theory', ...
       'NACA 0012 Vortex Panel', ...
       'NACA 2412 Vortex Panel', ...
       'NACA 4412 Vortex Panel', ...
       'Location', 'best');
grid on
xlim([-15 15])
set(gcf, 'Color', 'w')              % figure background white
set(gca, 'Color', 'w')              % axes background white

set(gca, 'XColor', 'k', 'YColor', 'k')   % axis lines + ticks black
set(gca, 'GridColor', 'k')               % grid lines black

set(get(gca,'Title'),'Color','k')        % title black
set(get(gca,'XLabel'),'Color','k')       % x-label black
set(get(gca,'YLabel'),'Color','k')       % y-label black
set(lgnd, 'Color', 'w', 'TextColor', 'k') % legend white with black text
hold off


%% Part 2
% This code uses Prandtl’s Lifting Line Theory (PLLT) to compute and plot
% the induced drag factor (delta) as a function of taper ratio for
% multiple aspect ratios, and compares results to Figure 5.20 from the
% Fundementals of Aerodynamics textbook.

%% Parameters
AR_vals = [4 6 8 10]; % Aspect ratios from textbook
lambda_vals = linspace(0.00001,1.0,100);% Taper ratio range

% Wing / airfoil properties (constant for validation)
c_r = 1;  % Root chord
a0_r = 2*pi; % Lift slope for root
a0_t = 2*pi; % Lift slope for tip
aero_r = 0; % Zero-lift angle for root
aero_t = 0; % Zero-lift angle for tip
geo_r = 5; % Geometric AoA for root
geo_t = 5; % Geometric AoA for tip
N = 50; % Number of Fourier terms

%% Plot  setup
figure(6);
hold on;
grid on;

%% Loop over aspect ratios
for k = 1:length(AR_vals)

    AR_target = AR_vals(k);
    delta_vals = zeros(size(lambda_vals));

    for i = 1:length(lambda_vals)

        lambda = lambda_vals(i);
        c_t = lambda * c_r;
        % Compute span for desired aspect ratio
        b_scaled = AR_target * (c_r + c_t) / 2;
        % Call PLLT function
        [e,~,~] = PLLT(b_scaled,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
        % Induced drag factor
        delta_vals(i) = 1/e - 1;

    end
    % Plot each AR curves on one plot
    plot(lambda_vals, delta_vals, 'LineWidth',2,'DisplayName',['AR = ' num2str(AR_target)])
end

%% Labels and formatting
xlabel('Taper ratio c_t / c_r')
ylabel('\delta')
title('Induced Drag Factor vs Taper Ratio (PLLT)')
legend('Location','northeast')
hold off;


%% Part 3 
%{
This part of the code uses the PLLF function to a wing of a Cessna 140. The
root airfoil is a NACA 2412 while the tip airfoil is a NACA 0012. 
%}

%% Parameters 
b3= 33+(4/12); % feet
c_r3= 5+(4/12); % feet
c_t3= 3+ (8.5/12); % feet
% Wing Area
S3= b3*(c_r3+c_t3)/2; % feet^2
AR3 = b3^2 / S3; 

%% NACA 0012 (tip) and  NACA 2412 (root) properties 

%Lift slope root and tip
% correct values, corrected from past submission
a0_r3= 6.88;
a0_t3= 6.88;
% Zero lift angle for root and tip
aero_t3 = 0; % degrees for symetrical airfoils
aero_r3= alphaL0_2412_vp; % degrees
alpha = 4; % degrees
%Geometric AoA for root and tip 
geo_r3= 1 + alpha; % degrees
geo_tVals= 0 + alpha; % degrees

%% Testing N values for convergence
N_vals3= 1:2:51; %odd terms (N = 51 gives up to four sig figs same values as N = 300)
CL_val3= zeros(size(N_vals3));
CDi_val3= zeros(size(N_vals3));

N_ref3= max(N_vals3); %reference the exact sln 

[e_3, CL_3,CDi_3]=PLLT(b3,a0_t3,a0_r3,c_t3,c_r3,aero_t3,aero_r3,geo_tVals,geo_r3,N_ref3);

delta_val= (1/e_3)-1;


 %% Computing CDi and Cl for each N value 
for i=1:length(N_vals3)

    N3=N_vals3(i);
   
    [e3,CL4,CDi3]= PLLT(b3,a0_t3,a0_r3,c_t3,c_r3,aero_t3,aero_r3,geo_tVals,geo_r3,N3);
    CL_val3(i)=CL4;
    CDi_val3(i)= CDi3;
end

error_Cl= abs(CL_val3-CL_3)./ CL_3*100;
error_CDi= abs(CDi_val3-CDi_3 )./ CDi_3*100;

target= [10 1 0.1];

for j=1:length(target)
total= target(j);

i_CL(j)=find( error_Cl < total,1,'first');
i_CDi(j)= find(error_CDi < total,1,'first');
end
%% Table for Deliverable 1
Percent_Error = ["10%"; "1%"; "0.1%"];
N_for_CL = [i_CL(1);  i_CL(2);  i_CL(3)];
N_for_CDi = [i_CDi(1); i_CDi(2); i_CDi(3)];
Values_for_CL = [CL_val3(i_CL(1));  CL_val3(i_CL(2));  CL_val3(i_CL(3))];
Values_for_CDi = [CDi_val3(i_CDi(1)); CDi_val3(i_CDi(2)); CDi_val3(i_CDi(3))];

Table1 = table(Percent_Error, N_for_CL, N_for_CL, Values_for_CL, Values_for_CDi);
disp(Table1)

%% Plots for Deliverable 2
%Cl convergence plot
figure(7);
plot(N_vals3,CL_val3,'-o');
hold on;

xline(N_vals3(find(error_Cl<10, 1)),'LineStyle','--', 'Color','r');
xline(N_vals3(find(error_Cl<1, 1)),'LineStyle','--','Color','g');
xline(N_vals3(find(error_Cl<0.1, 1)),'LineStyle','--','Color','b');
yline(CL_3, 'LineStyle','-');
legend('CL Values','10 Percent Error','1 Percent Error','0.1 Percent Error','Convergence Value');
xlabel('Odd terms');
ylabel('CL');
title('Convergence of CL');
grid on;

%CDi convergence plot
figure(8);
plot(N_vals3,CDi_val3,'-o');
hold on;

xline(N_vals3(find(error_CDi<10,1)),'LineStyle','--','Color','r');
xline(N_vals3(find(error_CDi<1,1)),'LineStyle','--','Color','g');
xline(N_vals3(find(error_CDi<0.1,1)),'LineStyle','--','Color','b');
yline(CDi_3,'LineStyle','-');
legend('CDi Values','10 Percent Error','1 Percent Error','0.1 Percent Error','Convergence Value');
xlabel('Odd terms');
ylabel('CDi');
title('Convergence of CDi');
grid on;

%% Deliverable 3
% using only 0.1 percent error so N = 8
% Lift force when V = 100 knots
V = 100 * 1.6878; % knots to ft/s conversion from unitconverters.net
% Standard Atmo numbers from Anderson Appendix E @ h = 10000ft
rho = 1.7556*10^(-3); %slugs / ft^3

Lift_Force = .5*rho*V^2*CL_val3(i_CL(3))*S3; % lbs
Induced_Drag_Force = .5*rho*V^2*CDi_val3(i_CDi(3))*S3; % lbs

% From Page 466 Theroy of Wing Sections @ Cl = 0.5284
cd = 0.007;

c_d_vals = [(0.008+0.0075)/2, (0.0075+0.007)/2, (0.0072+0.0061)/2,(0.007+0.006)/2, (0.0068+0.0059)/2, (0.0067+0.0059)/2, (0.0066+0.006)/2, (0.0067+0.0061)/2, (0.0068+0.0065)/2, (0.007+0.007)/2, (0.008+0.0069)/2];

% dyn_press = 0.5*rho*V^2;
% d_4 = c_d_lookup(end-1)*dyn_press*S3;
% Di_4 = c_Di_p(end)*dyn_press*S3;
% 
% L_4 = c_L_p(end)*dyn_press*S3;
% D_tot_4 = d_4 + Di_4;

Total_Drag_Force = .5*rho*V^2*(CDi_val3(i_CDi(3)) + cd)*S3;
Aerodynamic_Efficiency = Lift_Force / Total_Drag_Force;

Table2 = table(Lift_Force, Induced_Drag_Force, Aerodynamic_Efficiency);
disp(Table2)

%% Deliverable 4
% Total drag as a function of AoA
alpha_vals = linspace(-8,8,100); % range from -8 to 8 degrees

for i=1: length(alpha_vals)
    % only values that needed changed are geometric AoA
    geo_rVals= 1 + alpha_vals(i); % degrees
    geo_tVals= 0 + alpha_vals(i); % degrees
    [e_vals(i),CL_vals(i),CDi_vals(i)]= PLLT(b3,a0_t3,a0_r3,c_t3,c_r3,aero_t3,aero_r3,geo_tVals,geo_rVals,N3);
    Total_D_vals(i) = .5*rho*V^2*(CDi_vals(i) + cd)*S3;% lbs
    Lift_vals(i) = .5*rho*V^2*CL_vals(i)*S3; % lbs
    LOverD(i) = CL_vals(i) / (CDi_vals(i) + cd);
end
CD_total = CDi_vals + cd;
figure; 
plot(alpha_vals,CDi_vals,'b--');
hold on;
yline(cd,'r--');
hold on;
plot(alpha_vals,CD_total,'k');
xlabel('Angle of Attack (degrees)');
ylabel('Drag Coefficient');
title('Drag Coefficient vs Angle of Attack');
grid on;
legend('Induced Drag Coefficient','Sectional Drag Coefficient','Total Drag Coefficient','location','best');

%% Deliverable 5

figure; 
plot(alpha_vals,LOverD);
xlabel('Angle of Attack (degrees)');
ylabel('Lift Force / Total Drag Force');
title('Aerodynamic Efficiency of the Wing vs Angle of Attack');
grid on;


%% PLLT Function
function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%{
Inputs:
b = span
a0_t/r = lift slope (tip/root)
c_t/r = chord (tip/root)
aero_t/r = zero-lift angle (deg)
geo_t/r = geometric angle (deg)
N = number of Fourier terms

Outputs:
e = span efficiency factor
c_L = lift coefficient
c_Di = induced drag coefficient
%}

% Theta distribution
theta = zeros(1,N);
for i = 1:N
    theta(i) = i*pi/(2*N);
end

% Spanwise variation
c= c_r   + (c_t   - c_r  ).*cos(theta);
a0= a0_r  + (a0_t  - a0_r ).*cos(theta);
aero= aero_r+ (aero_t- aero_r).*cos(theta);
geo= geo_r + (geo_t - geo_r ).*cos(theta);

% Effective angle of attack
alpha_eff = deg2rad(geo - aero);

% Build system matrix
M = zeros(N,N);

for i = 1:N
    for j = 1:N
        n = 2*j - 1;
        M(i,j) = (4*b/(a0(i)*c(i))) * sin(n*theta(i)) ...
               + n*sin(n*theta(i))/sin(theta(i));
    end
end

% Solve for Fourier A coefficients
A = M \ alpha_eff';

% Compute delta
delta = 0;
for j = 2:N
    n = 2*j - 1;
    delta = delta + n*(A(j)/A(1))^2;
end

% Aerodynamic quantities
S  = b*(c_r + c_t)/2;
AR = b^2 / S;
e    = 1/(1 + delta);
c_L  = A(1)*pi*AR;
c_Di = c_L^2/(pi*e*AR);

end


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


%% Vortex Panel Function

function [CL] = Vortex_Panel_2(XB,YB,ALPHA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: %
% %
% XB = Boundary Points x-location %
% YB = Boundary Points y-location %
% VINF = Free-stream Flow Speed %
% ALPHA = AOA %
% %
% Output: %
% %
% CL = Sectional Lift Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
% Convert to Radians %
%%%%%%%%%%%%%%%%%%%%%%
ALPHA = ALPHA*pi/180;
%%%%%%%%%%%%%%%%%%%%%
% Compute the Chord %
%%%%%%%%%%%%%%%%%%%%%
CHORD = max(XB)-min(XB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the Number of Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = max(size(XB,1),size(XB,2))-1;
MP1 = M+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-Panel Relationships: %
% %
% Determine the Control Points, Panel Sizes, and Panel Angles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
IP1 = I+1;
X(I) = 0.5*(XB(I)+XB(IP1));
Y(I) = 0.5*(YB(I)+YB(IP1));
S(I) = sqrt( (XB(IP1)-XB(I))^2 +( YB(IP1)-YB(I))^2 );
THETA(I) = atan2( YB(IP1)-YB(I), XB(IP1)-XB(I) );
SINE(I) = sin( THETA(I) );
COSINE(I) = cos( THETA(I) );
RHS(I) = sin( THETA(I)-ALPHA );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships: %
% %
% Determine the Integrals between Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
for J = 1:M
if I == J
CN1(I,J) = -1.0;
CN2(I,J) = 1.0;
CT1(I,J) = 0.5*pi;
else
A = -(X(I)-XB(J))*COSINE(J) - (Y(I)-YB(J))*SINE(J);
B = (X(I)-XB(J))^2 + (Y(I)-YB(J))^2;
C = sin( THETA(I)-THETA(J) );
D = cos( THETA(I)-THETA(J) );
E = (X(I)-XB(J))*SINE(J) - (Y(I)-YB(J))*COSINE(J);
F = log( 1.0 + S(J)*(S(J)+2*A)/B );
G = atan2( E*S(J), B+A*S(J) );
P = (X(I)-XB(J)) * sin( THETA(I) - 2*THETA(J) ) ...
+ (Y(I)-YB(J)) * cos( THETA(I) - 2*THETA(J) );
Q = (X(I)-XB(J)) * cos( THETA(I) - 2*THETA(J) ) ...
- (Y(I)-YB(J)) * sin( THETA(I) - 2*THETA(J) );
CN2(I,J) = D + 0.5*Q*F/S(J) - (A*C+D*E)*G/S(J);
CN1(I,J) = 0.5*D*F + C*G - CN2(I,J);
CT2(I,J) = C + 0.5*P*F/S(J) + (A*D-C*E)*G/S(J);
CT1(I,J) = 0.5*C*F - D*G - CT2(I,J);
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships: %
% %
% Determine the Influence Coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
AN(I,1) = CN1(I,1);
AN(I,MP1) = CN2(I,M);
AT(I,1) = CT1(I,1);
AT(I,MP1) = CT2(I,M);
for J = 2:M
AN(I,J) = CN1(I,J) + CN2(I,J-1);
AT(I,J) = CT1(I,J) + CT2(I,J-1);
end
end
AN(MP1,1) = 1.0;
AN(MP1,MP1) = 1.0;
for J = 2:M
AN(MP1,J) = 0.0;
end
RHS(MP1) = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the gammas %
%%%%%%%%%%%%%%%%%%%%%%%%
GAMA = AN\RHS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Tangential Veloity and Coefficient of Pressure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
V(I) = cos( THETA(I)-ALPHA );
for J = 1:MP1
V(I) = V(I) + AT(I,J)*GAMA(J);
end
CP(I) = 1.0 - V(I)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Sectional Coefficient of Lift %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CIRCULATION = sum(S.*V);
CL = 2*CIRCULATION/CHORD;
end
