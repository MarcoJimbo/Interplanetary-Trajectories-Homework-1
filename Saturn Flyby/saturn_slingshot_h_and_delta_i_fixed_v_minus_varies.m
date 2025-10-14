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
% DeltaV = abs(V_helio_final - V_helio_planet)
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

r_p = 1.37e5;                       % 1) minimum flyby radius to not encounter the Saturn's rings[km]
v_minus = linspace(12,20,300);      % 2) heliocentric velocity magnitude of the SC @ SOI encounter [km/s]
delta_i = 45;                       % 3) Angle (in degrees) between SC velocity @SOI encounter and planet velocity [°]

% (NOTE: v_minus and delta_i depend on the interplanetary transfer arc from Earth to Saturn,
% they can either be assumed or evaluated in a separate script, together they form
% the VECTOR of the velocity of the SC @ encounter) 

%% Calculation of dependent flyby variables

v_inf = sqrt(v_minus.^2 + vp.^2 - 2*v_minus.*vp.*cosd(delta_i)); 
% Hyperbolic excess velocity (km/s), magnitude doesn't change, onl direction

S_D_i= (v_minus.*cosd(delta_i)-vp)/v_inf;
C_D_i= (v_minus.*sind(delta_i))/v_inf;
Delta_i = 2*atan(S_D_i./(1+C_D_i)); % initial angle between planet velocity and entrance hyperbolic excess velocity

e_hyp = 1 + (r_p/Mus)*v_inf.^2;   % eccentricity of the fly hyperbola
theta_star = acos(-1./e_hyp);     % asymptotic value of the true anomaly [°]
a_hyp = -Mus./(v_inf.^2);          % semi-major axis of the hyperbola [km]

xi = 2*theta_star - pi; % deflection angle: it's the angle that the VECTOR v_inf is rotated by [°]

% !!! WARNING !!!
% In a patched conics approximation there is ambiguity regarding the sense
% with which the flyby hyperbola is flown!!
% Therefore CHOOSE whether flyby is CW or CCW before proceeding !!!

% Final angle between planet velocity and exit hyperbolic excess velocity
  Delta_f = Delta_i + xi; % if CW
% Delta_f = Delta_i - xi; % if CCW

v_final = sqrt(v_inf.^2 + vp.^2 + 2*v_inf.*vp.*cos(Delta_f)); % velocity MAGNITUDE [km/s] in HC frame after flyby

C_d_f = (vp + v_inf.*cos(Delta_f))./v_final;
S_d_f = (v_inf.*sin(Delta_f))./v_final;
delta_f = 2*atan(S_d_f./(1+C_d_f));
Delta_V=abs(v_final-v_minus); %[km/s]
%% Plots and results

%% Table of the results
fprintf('Here below are reported the results in the case we fix the altitude of the spacecraft respect to the surface of the planet h and the angle between the velocity vecotr of the planet and the heliocentric velocity vector of the spacecraft at the encounter v-, delta_i \n while varying the heliocentric velocity vector of the spacecraft at the encounter\n\n');
i=1:length(v_minus);
v_minus_str=arrayfun(@(x) sprintf('%.2f km/s',x),v_minus,'UniformOutput',false);
v_inf_str=arrayfun(@(x) sprintf('%.2f km/s',x),v_inf,'UniformOutput',false);
xi_str=arrayfun(@(x) sprintf('%.2f°',x),rad2deg(xi),'UniformOutput',false);
v_final_str=arrayfun(@(x) sprintf('%.2f km/s',x),v_final,'UniformOutput',false);
delta_f_str=arrayfun(@(x) sprintf('%.2f°',x),rad2deg(delta_f),'UniformOutput',false);
Delta_V_str=arrayfun(@(x) sprintf('%.2f km/s',x),Delta_V,'UniformOutput',false);
T=table(i',v_minus_str',v_inf_str',xi_str',v_final_str',delta_f_str',Delta_V_str','VariableNames',{'#','Heliocentric magnitude velocity of SC','Initial hyperbolic excess velocity v_inf','deflection angle xi','Final heliocentric speed v_final','Final flight angle delta_f','Delta_V exiting Saturn'});
disp(T)


