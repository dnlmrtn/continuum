% [a b c] = Toroidal(0.2,2,2,1,1,1, 0)

a = 0.2;
r = 2;
l = 2;
g = 1;
mu = 1;
Q = 1;
% Define constants
A = 2*Q/(pi*(a^4));

% Define variables
syms rho phi s vx vy vz v_rho v_phi v_s temp x y z

% Define the fluid flow with respect to Toriodal coordinate system
v_rho = ((A^2)*g / (288*R*mu))*((a^2-rho^2)^2)*(4*(a^2)-(rho^2))*cos(phi);          % These equations are provided in the project description
v_phi = ((A^2)*g / (288*R*mu))*((a^2-rho^2)^2)*(4*(a^2)-23*(a^2)*(rho^2) + 7*(rho^4))*sin(phi);
v_s = A*((a^2)-(rho^2)) - ((3*A*rho*cos(phi))/(4*R))*((a^2) - (rho^2));
temp = diff(v_rho/sin(phi), rho);

p = (g*(A^2)*rho*((a^2)-(rho^2))*cos(phi)/R) + mu*cos(phi)*diff((temp + v_phi/(rho*sin(phi))-(v_rho/(rho*cos(phi)))),'rho',1);

% Convert to cartesian coordinate system (variables are still rho phi and s)
vx = (v_rho)*(cos(phi)*cos(s/R)) + (v_phi)*(-sin(phi)*cos(s/R)) + (v_s)*(-sin(s/R));        % These equations are provided in the project description
vy = v_rho*(sin(s/R)) + v_phi*sin(s/R) + v_s*(cos(s/R));
vz = v_rho*sin(phi) + v_phi*cos(phi);

% Define cartesian basis wrt toroidal coordinate system
x = (R + rho*cos(phi))*cos(s/R);           % These equations are provided in the project description
y = (R + rho*cos(phi))*sin(s/R);
z = rho*sin(phi);

%----- Plotting velocity field to verify results -----
% if plot == 1    
%     initalize some arrays
%     xSpace = [];
%     ySpace = [];
%     zSpace = [];
%     vxSpace = [];
%     vySpace = [];
%     vzSpace = [];
% 
%     for rho = 0:a/10:a
%         for phi = 0:(2*pi/10):(2*pi)
%             for s = 0:L/10:L
%                 x = (R + rho*cos(phi))*cos(s/R);
%                 y = (R + rho*cos(phi))*sin(s/R);
%                 z = rho*sin(phi);
%                 xSpace = [xSpace x];
%                 ySpace = [ySpace y];
%                 zSpace = [zSpace z];
%                 v_rho = ((A^2)*g / (288*R*mu))*((a^2-rho^2)^2)*(4*(a^2)-(rho^2))*cos(phi);
%                 v_phi = ((A^2)*g / (288*R*mu)) * ((a^2-rho^2)^2) * (4*(a^2)-23*(a^2)*(rho^2) + 7*(rho^4))*sin(phi);
%                 v_s = A*(a^2-rho^2) - ((3*A*rho*cos(phi))/(4*R))*(a^2 - rho^2);
%                 vx = (v_rho)*(cos(phi)*cos(s/R)) + (v_phi)*(-sin(phi)*cos(s/R)) + (v_s)*(-sin(s/R));
%                 vy = v_rho*(sin(s/R)) + v_phi*sin(s/R) + v_s*(cos(s/R));
%                 vz = v_rho*sin(phi) + v_phi*cos(phi);
%                 vxSpace = [vxSpace vx];
%                 vySpace = [vySpace vy];
%                 vzSpace = [vzSpace vz];
%             end
%         end
%     end
% end
% quiver3(xSpace,ySpace,zSpace,vxSpace,vySpace,vzSpace)

% Now that we have cartesian velocity field in terms of toroidal
% coordinates, we can do integration


% Rate of deformation tensor:
gradV = jacobian([vx vy vz], [rho phi s])*(jacobian([x,y,z],[rho,phi,s]).^-1);
D = 0.5*(gradV+gradV');
T = -p*eye(3) + 2*mu*D;

%Define normal of pipe surface wrt cartesian coordinate system
syms n
n = [(rho*cos(phi))*cos(s/R);(rho*cos(phi))*sin(s/R); rho*sin(phi)];

% 2D surface integral over the surface of the pipe
rho = subs(rho,a);
[Fx, Fy, Fz] = int(int(T*n,phi,0,2*pi),s,0,L);