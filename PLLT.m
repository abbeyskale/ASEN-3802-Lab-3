function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%{
Inputs: 
b is span (in feet)
a0_t is the cross-sectional lift slope at the tips (per radian)
a0_r is the cross-sectional lift slope at the root (per radian)
c_t is the chord at the tips (in feet)
c_r is the chord at the root (in feet)
aero_t is the zero-lift angle of attack at the tips (in degrees)
aero_r is the zero-lift angle of attack at the root (in degrees)
geo_t is the geometric angle of attack at the tips (in degrees)
geo_r is the geometric angle of attack at the root (in degrees)
N is number of terms in series expansion

Outputs: 
e = space efficiancy factor
C_L = coeff of lift
C_Di = induced coeff of drag
%}

% solve linearly
y = linspace(0,b/2,100);
for i= 1:size(y)
if y(i) > -b/2 && y(i) < 0
    c = c_r + (c_r-c_t)*y/(b/2);
    a0 = a0_r + (a0_r-a0_t)*y/(b/2);
    aero = aero_r + (aero_r-aero_t)*y/(b/2);
    geo = geo_r + (geo_r-geo_t)*y/(b/2);
else 
    c = c_r + (c_t-c_r)*y/(b/2);
    a0 = a0_r + (a0_t-a0_r)*y/(b/2);
    aero = aero_r + (aero_t-aero_r)*y/(b/2);
    geo = geo_r + (geo_t-geo_r)*y/(b/2);
end
end

% left hand side is effective angle of attack
alpha_eff = geo - aero;
n_odd = 1:2:N;

M = zeros(N,N);
theta = zeros(1,N);
for i=1:N
    theta(i) = i*pi/(2*N);
end
for i=1:size(a0)
    for j=1:size(a0)
    n = n_odd(j);
    test(i,j) = sin(n*theta(i));
    test2(i,j) = 4*b/(a0(i)*c(i));
    test3(i,j) = n/sin(theta(i));
    M(i,j) = sin(n*theta(i)) * (4*b./(a0(i)*c(i)) + n/sin(theta(i)));
    end
end

A = M ./ alpha_eff;

delta = 0;

% for j = 2:N
%     n = n_odd(j);
%     delta = delta + n*(A(i)/A(1))^2;
% end
% 
% e = 1 / (1 + delta);
% 
% S = b*((c_r - c_t)/2 + c_t);
% c_L = A(1) * pi * b^2 / S;

