%% Uranus arrival
% Once we have exited the Saturn SOI we need to arrive to the Uranus SOI so
% we have to analyze the transfer arc from Saturn to Uranus.
clc 
close all

%% Definition of the problem constants (all in km and km/s)
MuS = 1.32712440018e11;         % Sun GM [km^3/s^2]
MuU = 5793939;                  % Uranus GM [km^3/s^2]
Ru = 2872.46e6;                 % Uranus semimajore axis [km]
Rs  = 9.537e8;                  % Saturn semimajor axis [km]
r_soi = Ru*(MuU/MuS)^(2/5);     % Uranus Sphere of Influence radius [km]

%% Transfer arc variable to go from Saturn to Uranus (1 independent)
v_minus = 18;                       % S/C heliocentric velocity outside the Saturn SOI [km/s]
a = (2/Rs - v_minus^2 / MuS)^-1;    % semimajor axis of interplanetary transfer arc [km]
e = 1-Rs/a;                         % eccentricity of interplanetary transfer arc

%% True anomaly and velocity at the arrival of Uranus SOI
value=(a*(1-e^2)/(e*Ru))- 1/e;                                  % variable to check the value of true anomaly
theta_star_plus = acos(value);                                  % true anomaly of S/C outside the Uranus SOI
v_r_plus = sqrt(MuS/a*(1-e^2))*e*sin(theta_star_plus);          % radial orbital component of S/C velocity at the encounter of Uranus SOI
v_theta_plus = sqrt(MuS/a*(1-e^2)*(1+e*cos(theta_star_plus)));  % tangent orbital component of S/C velocity at the encounter of Uranus SOI
v_plus = sqrt(v_r_plus^2 + v_theta_plus^2);

%% Check if the value of true anomaly is in the range [-1,1]

if value<-1 || value>1
    fprintf('ERROR!! the true anomaly is a complex number, check the heliocentric S/C velocity');
    return
end
%% Check if the eccentricity of the transfer arc is in the range [0,inf]

if e<0
    fprintf('ERROR!! the eccentrcity must be a positive number, check the heliocentric S/C velocity');
    return
end
%% Plots and results

fprintf('True anomaly at the encounter of Uranus SOI : %.2f°\n', rad2deg(theta_star_plus));
fprintf('The S/C heliocentric velocity at the encounter of Uranus SOI : %.2f km/s\n', v_plus);

%% Plot of the trasfer arc in the Sun centered frame
figure; hold on; axis equal; grid on
title('Uranus Transfer Arc (Sun-centered frame)');
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

% Draw Uranus orbit for context
theta_circle = linspace(0,2*pi,200);
plot(Ru*cos(theta_circle), Ru*sin(theta_circle), '-c');

legend('Interplanetary transfer arc','Sun','Saturn orbit','Uranus orbit');