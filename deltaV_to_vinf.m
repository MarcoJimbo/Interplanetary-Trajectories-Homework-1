function v_inf = deltaV_to_vinf(delta_V,h_leo)
% Function to evaluate the hyperbolic excess velocity (v_inf) from an impulsive velocity change (delta_V) applied in a Low Earth Orbit (LEO) at altitude h_LEO.
%
% INPUT: delta_V [km/s] and h_LEO [km]
%
% OUTPUT: v_inf [km/s]  
%
%% DATA
mu_E = 398600; % [km^3/s^2] 
mu_sun= 1.327*10^11; % [km^3/s^2] 
r_e_sun= 149.6*10^6;
R_E  = 6378;         % [km] Earth radius         
r_leo = R_E + h_leo; % [km]
%%
% Circular orbital velocity in LEO
v_c = sqrt(mu_E / r_leo);

% Velocity after the applied delta-V
v_p = v_c + delta_V;

% Escape velocity from that altitude
 v_esc = sqrt(2 * mu_E / r_leo);


% Check condition
if v_p < v_esc
        error('Delta-V too small: trajectory is not hyperbolic (v_p < v_esc)');
end

% Calculation of hyperbolic excess velocity (v_inf)
v_inf = sqrt(v_p^2 - 2 * mu_E / r_leo);

e_hyp=1+(r_leo/mu_E)*(v_inf)^2;

v_E=sqrt(mu_sun/r_e_sun);

end
