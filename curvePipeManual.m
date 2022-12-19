% Define v wrt toroidal coordinate system
a = 1;
R = 1;
L = 1;
g = 2;
mu = 1;
Q = 1;
rho = a;

A = 2*Q/(pi*(a^4));
syms phi s 
% e_rho e_phi e_s v_rho v_phi v_s temp Fx Fy Fz

% Velocity field
v_rho = ((A^2)*g / (288*R*mu)) * ((a^2-rho^2)^2)*(4*(a^2)-(rho^2))*cos(phi);
v_phi = ((A^2)*g / (288*R*mu)) * ((a^2-rho^2)^2) * (4*(a^2)-23*(a^2)*(rho^2) + 7*(rho^4))*sin(phi);
v_s = A*(a^2-rho^2) - ((3*A*rho*cos(phi))/(4*R))*(a^2 - rho^2);

% Hand calculated partial derivatives
dvrho_drho = -(A^2)*g*cos(phi)*rho*(rho^2-3*(a^2))*((rho^2)-(a^2))/(48*R*mu);
dvrho_dphi = (A^2)*g*((rho^2)-4*(a^2))*((rho^2 - a^2)^2)*sin(phi);
dvrho_ds = 0;

dvphi_drho = -(A^2)*g*sin(phi)*rho*(21*rho^4 - 60*(a^2)*(rho^2) + 23*a^4 + 4*a^4)/(144*R*mu);
dvphi_dphi = (A^2)*g*(a^2 - rho^2)*(7*rho^4 - 23*(a^2)*(rho^2) + 4*a^2)*cos(phi)/(288*R*mu);
dvphi_ds = 0;

dvs_drho = A*(9*cos(phi)*(rho^2) - 8*R*rho - 3*(a^2)*cos(phi))/(4*R);
dvs_dphi = 3*A*rho*(a^2 - rho^2)*sin(phi)/(4*R);
dvs_ds = 0;

% Pressure field
p = (g*(A^2)*rho*((a^2)-(rho^2))*cos(phi)/R) + mu*cos(phi)*diff((dvrho_drho/sin(phi) + v_phi/(rho*sin(phi))-(v_rho/(rho*cos(phi)))),'rho',1);

% Gradient of velocity field
gradV = [dvrho_drho dvrho_dphi dvrho_ds;
         dvphi_drho dvphi_dphi dvphi_ds;
         dvs_drho dvs_dphi dvs_ds];

% Calculate T
D = 0.5*(gradV + gradV.');
T = -p*eye(3) + 2*mu*D;
n = [1; 0; 0];

% Calculate F
F = T * n;

Fx = F(1)*cos(phi)*cos(s/R)*F(1) - F(2)*sin(phi)*cos(s/R) - F(3)*sin(s/R);
Fy = F(1)*cos(phi*sin(s/R)) - F(2)*sin(phi)*cos(s/R) + F(3)*cos(s/R);
Fz = F(1)*sin(phi) + F(2)*cos(phi);

force_cart = [Fx Fy Fz];
%force_cart = subs(force_cart, rho, a);

ringForce = int(force_cart,phi, 0, 2*pi);

totalForce = int(ringForce,s,0,L);

double(totalForce)


