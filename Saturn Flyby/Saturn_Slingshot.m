%% SATURN SLINGSHOT ANALYSIS 

% As part of a interplanetary mission towards Uranus, the option of 
% exploiting a gravitational assist around Saturn is thereby evalued from
% an energetic POV. 

%% Script Description:
% This script calculates the velocity magnitude & vector of
% the spacecraft (in the Heliocentric frame) after a flyby of Saturn
% for a wide range of flyby variables (the periapsis radius r_p,
% and the spacecraft SOI entrance velocity)
%
% With this the Delta V of the gravitational slingshot is evalued 
%
%% Theoretical background:
% In a patched-conics a flyby trajectory is a hyperbola.
% The trajectory rotates the so called hyperbolic excess velocity vector by angle ξ (xi), 
% thus producing a velocity change in the Heliocentric reference frame,
% indeed: 
% V_helio_final = V_helio_planet + V_inf_final ==> 
%
% DeltaV = abs(V_helio_final - V_helio_planet) = abs(V_inf_final)
%
% This DeltaV constitutes the gravitational slingshot effect

clc
close all
%% Definition of the problem constants (all in km and km/s)
Mus = 3.7931187e7;       % Saturn GM [km^3/s^2]
MuS = 1.32712440018e11;  % Sun GM [km^3/s^2]
rs  = 58232;             % Saturn equatorial radius [km]
Rs  = 9.537e8;           % Saturn semimajor axis ≈ 9.537 AU [km]
vp  = sqrt(MuS/Rs);      % Saturn orbital speed [km/s] ≈ 9.68 km/s
r_soi = Rs*(Mus/MuS)^(2/5); % Sphere of Influence radius [km]


%% Problem INPUTS (aka the flyby variables, 3 independent ones)
h = 30;                  % flyby altitude [km]
r_p = h + rs;            % 1) flyby radius [km]
v_minus = 12;            % 2) heliocentric velocity magnitude of the SC @ SOI encounter [km/s]
delta_i = 30;            % 3) Angle (in degrees) between SC velocity @SOI encounter and planet velocity

% (NOTE: v_minus and delta_i depend on the interplanetary transfer arc from Earth to Saturn,
% they can either be assumed or evaluated in a separate script, together they form
% the VECTOR of the velocity of the SC @ encounter) 

%% Calculation of dependent flyby variables

v_inf = sqrt(v_minus^2 + vp^2 - 2*v_minus*vp*cosd(delta_i)); 
% Hyperbolic excess velocity (km/s), magnitude doesn't change, only direction

S_D_i= (v_minus*cosd(delta_i)-vp)/v_inf;
C_D_i= (v_minus*sind(delta_i))/v_inf;
Delta_i = 2*atan(S_D_i/(1+C_D_i)); % initial angle between planet velocity and entrance hyperbolic excess velocity

e_hyp = 1 + (r_p/Mus)*v_inf^2;   % eccentricity of the fly hyperbola
theta_star = acos(-1/e_hyp);     % asymptotic value of the true anomaly
a_hyp = -Mus/(v_inf^2);          % semi-major axis of the hyperbola [km]

xi = 2*theta_star - pi; % deflection angle: it's the angle that the VECTOR v_inf is rotated by 

% !!! WARNING !!!
% In a patched conics approximation there is ambiguity regarding the sense
% with which the flyby hyperbola is flown!!
% Therefore CHOOSE whether flyby is CW or CCW before proceeding !!!

% Final angle between planet velocity and exit hyperbolic excess velocity
  Delta_f = Delta_i + xi; % if CW
% Delta_f = Delta_i - xi; % if CCW

v_final = sqrt(v_inf^2 + vp^2 + 2*v_inf*vp*cos(Delta_f)); % velocity MAGNITUDE [km/s] in HC frame after flyby

C_d_f = (vp + v_inf*cos(Delta_f))/v_final;
S_d_f = (v_inf*sin(Delta_f))/v_final;
delta_f = 2*atan(S_d_f/(1+C_d_f));

%% Plots and results

fprintf('Initial hyperbolic excess velocity: %.2f km/s\n', v_inf);
fprintf('Deflection angle xi: %.2f deg\n', rad2deg(xi));
fprintf('Final heliocentric speed: %.2f km/s\n', v_final);
fprintf('Final flight angle delta_f: %.2f deg\n', rad2deg(delta_f));

%% Plot of the flyby trajectory in Saturn-centered frame

figure; hold on; axis equal; grid on
title('Saturn Flyby Trajectory (Saturn-centered frame)');
xlabel('x [km]'); ylabel('y [km]');

% Range of true anomaly
theta = linspace(-theta_star, theta_star, 1000);
r = (a_hyp*(e_hyp^2 - 1))./(1 + e_hyp*cos(theta));

% Convert to Cartesian
x = r .* cos(theta);
y = r .* sin(theta);

% Plot hyperbolic trajectory
plot(x, y, 'b', 'LineWidth', 1.5);

% Plot periapsis point (+rp or -rp if cw o ccw)
plot(-r_p, 0, 'ro', 'MarkerFaceColor','r');

% Draw Saturn with PNG texture
[img,~,alpha] = imread('saturn.png');
saturn_radius_plot = rs; % [km]
xRange = [-1 1]*saturn_radius_plot;
yRange = [-1 1]*saturn_radius_plot;
image('XData',xRange,'YData',yRange,...
      'CData',img,'AlphaData',alpha);

% Draw SOI for context
theta_circle = linspace(0,2*pi,200);
plot(r_soi*cos(theta_circle), r_soi*sin(theta_circle), '--k');

legend('Hyperbolic trajectory','Periapsis','Saturn','SOI');
