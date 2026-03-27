function [x_b,y_b,x_cam, y_cam] = NACA_Airfoils(m,p,t,c,N)

%Thicknesss distribution 
yt= (c*t/0.2)*((0.2969*sqrt(x/c)) - (0.126*(x/c))-(0.3516*(x/c)^2)+(0.2843*(x/c)^3)-(0.1036*(x/c)^4));

% camber lineand slope 
if p==0
yc= zeros(size(x));
dyc_dx= zeros(size(x));
end

% For NACA 2421
%if(x>=0) && (x<pc)
   % yc= m*(x/p^2)*(2*p-(x/c));
%elseif (x>= pc) && (c<=x)
% yc= m*(((c-x)/(1-p)^2)*(1+(x/c)-2*p));
%else
%    fprintf('Error');
%end

if 
end

% Angle
Xi= arctan(dyc_dx);

% Upper Surface
xU= x-yt*sin(Xi);
yU= yc+yt*cos(Xi);

%Lower Surface
xL= x+yt*sin(Xi);
yL= yc+yt*cos(Xi);
end
