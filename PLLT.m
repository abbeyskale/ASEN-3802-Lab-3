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
for i=1:N
    theta(i) = i*pi/(2*N);
end

c = c_r + (c_t-c_r)*cos(theta);
a0 = a0_r + (a0_t-a0_r)*cos(theta);
aero = aero_r + (aero_t-aero_r)*cos(theta);
geo = geo_r + (geo_t-geo_r)*cos(theta);

% left hand side is effective angle of attack
alpha_eff = geo - aero;

M = zeros(N,N);
theta = zeros(1,N);

for i=1:length(N)
    for j=1:length(N)
    n = 2*j - 1;
    M(i,j) = (4*b/(a0(i)*c(i)) * sin(n*theta(i))) + n*sin(n*theta(i))/sin(theta(i));
    end
end

A = M \ alpha_eff';

delta = 0;
e = 0;
c_L = 0;
c_Di = 0;

for j = 2:N
    n = 2*j -1;
    delta = delta + n*(A(i)/A(1))^2;
end

S = b*((c_r - c_t)/2 + c_t);
AR = b^2 / S;
e = 1 / (1 + delta);
c_L = A(1)*pi*AR;
c_Di = c_L^2 / (pi*e*AR);
