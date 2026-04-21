clc; clear; close all;
%% Part 3

%% Parameters
% Wing / airfoil properties (constant for validation)
c_r = 5.444;  % Root chord in inches
a0_r = 2*pi; % Lift slope for root
a0_t = 2*pi; % Lift slope for tip
aero_r = 5; % Zero-lift angle for root
aero_t = 5; % Zero-lift angle for tip
geo_r = 1; % Geometric AoA for root in degrees
geo_t = 0; % Geometric AoA for tip in degrees
N = 10000; % Number of Fourier terms

b = 33.333; % feet
c_t = 3.7083; % feet

S  = b*(c_r + c_t)/2;
AR = b^2 / S;

[e,C_L,C_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,...
                aero_t,aero_r,geo_t,geo_r,N);
% Induced drag factor
delta_val = 1/e - 1;

%% PLLT Function
function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,...
                            aero_t,aero_r,geo_t,geo_r,N)

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
c    = c_r   + (c_t   - c_r  ).*cos(theta);
a0   = a0_r  + (a0_t  - a0_r ).*cos(theta);
aero = aero_r+ (aero_t- aero_r).*cos(theta);
geo  = geo_r + (geo_t - geo_r ).*cos(theta);

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
