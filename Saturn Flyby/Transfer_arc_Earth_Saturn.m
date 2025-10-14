%% Earth departure analysis
% To analyze the interplanetary mission toward Uranus moons we have to start to analyze the initial step. 
% We want to know the parameters that characterize the interplanetary transfer arc from Earth to Saturn considering
% the spacecraft outside the SOI of Earth with a heliocentric velocity v_plus which is the sum of the hyperbolic excess velocity v_inf
% and the heliocentric orbital velocity of Earth ve.
clc
close all

%% Definition of the problem constants (all in km and km/s)
MuS = 1.32712440018e11;         % Sun GM [km^3/s^2]
MuE = 398.6;                    % Earth GM [km^3/s^2]
Re = 1.49380e8;                 % Earth semimajore axis [km]
Rs  = 9.537e8;                  % Saturn semimajor axis [km]
ve = 29.8;                      % Earth orbital speed [km/s]

%% Transfer arc variable to go from Earth to Saturn (1 independent)
% When we are outside the Earth SOI the hyperbolic excess velocity v_inf
% coincides with the difference between the S/C velocity at perigee of the 
% Earth hyperbola v_plus, and the heliocentric orbital velocity of Earth
% ve. So in this case v_inf coincides with the energetic cost to escape
% Earth and do a transfer arc without power.

v_inf = 14;                                       % hyperbolic excess velocity [km/s] (minimum and maximum velocity to have an ellipse : 9.399 12.352)

%% Parameters that characterize the transfer arc near the Earth
v_minus = v_inf + ve;                                 % heliocentric velocity outside the Erath SOI [km/s]
a = (2/Re - v_minus^2 / MuS)^-1;                       % semimajor axis of interplanetary transfer arc [km]
e = 1-Re/a;                                          % eccentricity of interplanetary transfer arc
%% True anomaly and velocity at the arrival of Saturn SOI
value=(a*(1-e^2)/(e*Rs))- 1/e;                                  % variable to check the value of true anomaly
theta_star_plus = acos(value);                                  % true anomaly of S/C outside the Saturn SOI
v_r_plus = sqrt(MuS/a*(1-e^2))*e*sin(theta_star_plus);          % radial orbital component of S/C velocity at the encounter of Saturn SOI
v_theta_plus = sqrt(MuS/a*(1-e^2)*(1+e*cos(theta_star_plus)));  % tangent orbital component of S/C velocity at the encounter of Saturn SOI
v_plus = sqrt(v_r_plus^2 + v_theta_plus^2);

%% Check if the value of true anomaly is in the range [-1,1]

if value<-1 || value>1
    fprintf('ERROR!! the true anomaly is a complex number, check the hyperbolic excess velocity');
    return
end
%% Check if the eccentricity of the transfer arc is in the range [0,inf]

if e<0
    fprintf('ERROR!! the eccentrcity must be a positive number, check the hyperbolic excess velocity');
    return
end
%% Plots and results

fprintf('True anomaly at the encounter of Saturn SOI : %.2f°\n', rad2deg(theta_star_plus));

%% Plot of the trasfer arc in the Sun centered frame
figure; hold on; axis equal; grid on
title('Saturn Transfer Arc (Sun-centered frame)');
xlabel('x [km]'); ylabel('y [km]');

% Range of true anomaly to plot the transfer arc
theta = linspace(0, theta_star_plus, 1000);
r = (a*(1-e^2))./(1 + e*cos(theta));

% Convert to Cartesian
x = r .* cos(theta);
y = r .* sin(theta);

% Plot transfer arc trajectory
plot(x, y, 'b', 'LineWidth', 1.5);

% Plot Sun point
plot(0, 0, 'ro', 'MarkerFaceColor','y');

% Draw Saturn orbit for context
theta_circle = linspace(0,2*pi,200);
plot(Rs*cos(theta_circle), Rs*sin(theta_circle), '-m');

% Draw Earth orbit for context
theta_circle = linspace(0,2*pi,200);
plot(Re*cos(theta_circle), Re*sin(theta_circle), '-g');


legend('Interplanetary transfer arc','Sun','Saturn orbit','Earth orbit');
