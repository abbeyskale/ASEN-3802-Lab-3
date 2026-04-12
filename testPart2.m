clc; clear; 

b = 40; %feet
a0_t = 1;
a0_r = 1.5;
c_t = 5; 
c_r = 20;
aero_t = 10;
aero_r = 10;
geo_t = 10;
geo_r = 10;
N = 3;

[e,cl,cdi] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)