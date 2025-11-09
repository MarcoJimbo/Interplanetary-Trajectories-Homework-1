%% User-defined constants and functions
% NY = number of differential equations
% KP = number of unknown parameters (not initial conditions)
% K = = total number of unknowns (parameters + initial conditions)

% params = structure containing all problem data that are used in the
% USER_DEFINED function

%% USER_DEFINE functions
% The code requires you to define the functions
% - setIC
% - funz
% - bound 
% - discon
% - post_process
% 
% function [val, isExplicit] = setIC(j)
% Set explicit value of some initial conidition s(j)
% if s(j) == val, isExplicit must be true
% 




function randevouz_tfNOTfree
clear 
close all
clc

prb.PBIS = 1;      % Bisection coefficient (>1) 
prb.RMIN = 1e-1;    % Coefficient of Newton method: x_{k+1} = x_{x} - RMIN * f(x_k)/Df_{k} 
prb.JMAX = 100;    % Maximum number of steps 

prb.NY = 10;       % Number of differential equations
prb.KP = 3;        % Number of unknowns that are not initial boundary values (e.g., the time-length of an arc)
prb.K = 5+3;       % Number of unknowns (total). Must be equal to the number of boundary constraints
prb.NCMP = 3;      % Number of arcs (>=1) 
prb.ATOL = 1e-6;   % Absolute tolerance for ODE solver
prb.RTOL = 1e-8;   % Relative tolerance for ODE solver

params.mu = 1;                    % Gravitational Parameters
params.Q = 0.2;                   % Max flux
params.C = 3;                     % Exhaust velocity
params.THR = params.Q*params.C;   % Thrust

prb.params = params;
prb.setIC = @setIC;
prb.funz = @funz;
prb.discon = @discon;
prb.bound = @bound;

prb.fsolve_opz = optimset('Display','Iter');
% fsolve_opz.Algorithm = 'levenberg-marquardt';
prb.fsolve_opz.MaxIter = 100; 
prb.fsolve_opz.MaxFunEvals=1e6;
prb.solve = 1; 

% TOF = 3.5 
%[tau1,tau2,tau3,λr(0),λtheta(0),λvr(0),λvt(0),λm(0)]
% prb.YP_guess = [0.5 5 0.5 0.5 -0.5 -0.5 0.5 0.5]';
% prb.YP_guess = [0.4 2.6 0.25 -0.35 -0.4 -0.6 0.5 0.03]';
% prb.YP_guess = [0.3578    2.6469    0.2692    0.0986   -0.5131   -0.4647    0.6504   -0.5629]';
% prb.YP_guess = [0.470674  2.604793  0.340936  0.748405  -0.301882  0.658846  0.822139  -1.221989]';
% prb.YP_guess = [0.628962  2.447314  0.392816  0.750363  -0.113146  0.587829  0.688895  -1.033808]';
% prb.YP_guess = [0.705697  2.349251  0.440869  0.367790  -0.036985  0.303248  0.388820  -0.715674]';
% prb.YP_guess = [0.714547  2.336518  0.448000  0.282819  -0.028495  0.240299  0.315521  -0.626729]';
% prb.YP_guess = [0.716971  2.332883  0.450129  0.255277  -0.025962  0.219730  0.291047  -0.595268]';
% prb.YP_guess = [0.717015  2.332816  0.450169  0.254752  -0.025914  0.219337  0.290578  -0.594655]';
prb.YP_guess = [1.072358  1.002614  0.825045  0.596054  -0.013899  0.363573  0.446344  -0.224054]';




sol = solve_by_shooting(prb);
prb.YP_guess=sol.YP;

post_process(prb,sol, params)

end

 


%% Problem dependant functions
% S = [r, theta, vr, vt, m, λr, λtheta, λvr, λvt, λm]
%     [1,     2,  3,  4, 5,  6,      7,   8,   9, 10]

function [val, isExplicit] = setIC(j)
% Set explicit value of some initial conidition s(j)
% if s(j) == val, isExplicit must be true

switch j
case 1 % r, 
    isExplicit = 1; val = 1;
case 2 % theta
    isExplicit = 1; val = 0;    
case 3 % vr, 
    isExplicit = 1; val = 0;
case 4 % vt
    isExplicit = 1; val = 1;
case 5 % m
    isExplicit = 1; val = 1; 
otherwise
    val = NaN;
    isExplicit = false;
end
end

function ds = funz(t, s, icmp, YP, params)
% ODE function to be integrated.
% S = [r, theta, vr, vt, m, λr, λtheta, λvr, λvt, λm]
%     [1,     2,  3,  4, 5,  6,      7,   8,   9, 10]

THR = params.THR;
C = params.C;
mu = params.mu;

tau = YP(icmp);

% extract from state
r = s(1);
theta = s(2);
vr = s(3);
vt = s(4);
m = s(5);

lambda_r = s(6);
lambda_theta = s(7);
lambda_vr = s(8);
lambda_vt = s(9);
lambda_m = s(10);


% % Check
% if any(~isfinite(s)) || r < 0 || m < 0
%     error('Variabili non fisiche: r=%.2f, m=%.2f', r, m)
% end

% Optimal control
if icmp==1 || icmp==3
    phi = atan2(-lambda_vt, -lambda_vr);
    % phi = atan2(lambda_vt, lambda_vr);
    Tr = THR*cos(phi);
    Tt = THR*sin(phi);
    dmdt = - THR / C;
else % Coast (icmp == 2)
    Tr = 0; 
    Tt = 0; 
    dmdt = 0;
end

% if (m<0.05)
%     Tr = 0; 
%     Tt = 0;
%     dmdt=0;
% end

    % state dynamics
    drdt = vr;
    dthetadt = vt/r;
    dvrdt = (vt^2)/r - mu/r^2 + Tr/m;
    dvtdt = (-vr*vt)/r + Tt/m;


    % costates dynamics
    dlambda_r = lambda_vr*(vt^2 - 2*mu/r)/r^2 - lambda_vt*(vr*vt)/r^2 + lambda_theta*(vt/r^2);
    dlambda_theta = 0;
    dlambda_vr = -lambda_r + lambda_vt*vt/r;
    dlambda_vt = -2*lambda_vr*vt/r + lambda_vt*vr/r - lambda_theta/r;
    dlambda_m = (lambda_vr*Tr + lambda_vt*Tt)/m^2;
    
    ds = tau*[drdt; dthetadt; dvrdt; dvtdt; dmdt; dlambda_r; dlambda_theta; dlambda_vr; dlambda_vt; dlambda_m];
end


function swf = switchingF(t,s,icmp,YP, params)
% Given the state, compute the Hamiltonian function
% S = [r, theta, vr, vt, m, λr, λtheta, λvr, λvt, λm]
%     [1,     2,  3,  4, 5,  6,      7,   8,   9, 10]

C=params.C;

% extract from state
m = s(5);

lambda_vr = s(8);
lambda_vt = s(9);
lambda_m = s(10);

swf = -sqrt(lambda_vr^2 + lambda_vt^2)/m - lambda_m/C;

end



function H = hamiltonian(t,s,icmp,YP, params)
% Given the state, compute the Hamiltonian function
% S = [r, theta, vr, vt, m, λr, λtheta, λvr, λvt, λm]
%     [1,     2,  3,  4, 5,  6,      7,   8,   9, 10]
mu=params.mu;
THR=params.THR;

% extract from state
r = s(1);
theta = s(2);
vr = s(3);
vt = s(4);
m = s(5);

lambda_r = s(6);
lambda_theta = s(7);
lambda_vr = s(8);
lambda_vt = s(9);
lambda_m = s(10);

%optimal control
if icmp==1 || icmp==3
    phi = atan2(-lambda_vt, -lambda_vr);
    % phi = atan2(lambda_vt, lambda_vr);
    Tr = THR*cos(phi);
    Tt = THR*sin(phi);
    dmdt = -THR / params.C;
else % Coast (icmp == 2)
    Tr = 0; 
    Tt = 0; 
    dmdt = 0;
end


% f vector
drdt = vr;
dthetadt = vt/r;
dvrdt = (vt^2)/r - mu/r^2 + Tr/m;
dvtdt = (-vr*vt)/r + Tt/m;
f = [drdt;dthetadt;dvrdt;dvtdt;dmdt];

% Lambda vector
lambda = [lambda_r;lambda_theta;lambda_vr;lambda_vt;lambda_m];

% Hamiltonian (L=0)
H = lambda'*f;

end


function ER = bound(S_left, S_right, YP, params)
% Computes the error between the final values (from integration)
% and the final boundary condition.
% S = [r, theta, vr, vt, m, λr, λtheta, λvr, λvt, λm]
%     [1,     2,  3,  4, 5,  6,      7,   8,   9, 10]

    % Condizioni finali (t=7)
    s_final = S_right(:, end); 
    r_final = s_final(1);
    theta_final = s_final(2);
    vr_final = s_final(3);
    vt_final = s_final(4);
    lambda_r_final = s_final(6);
    lambda_theta_final = s_final(7);
    lambda_vr_final = s_final(8);
    lambda_vt_final = s_final(9);
    lambda_m_final = s_final(10);
    lambda_norm_sq = lambda_r_final^2+lambda_theta_final^2+lambda_vr_final^2+lambda_vt_final^2+lambda_m_final^2;

    swf1 = switchingF(1, S_right(:,1), 1, YP, params);
    swf2 = switchingF(1, S_right(:,2), 2, YP, params); 
    swf3 = switchingF(1, S_right(:,3), 3, YP, params); 
    H1minus = hamiltonian(1,S_right(:,1),1,YP,params);
    H2minus = hamiltonian(1,S_right(:,2),2,YP,params);
    H12 = hamiltonian(1,S_right(:,1),1,YP,params) - hamiltonian(1,S_left(:,2),2,YP,params);
    H23 = hamiltonian(1,S_right(:,2),2,YP,params) - hamiltonian(1,S_left(:,3),3,YP,params);
    H3final = hamiltonian(1,S_right(:,3),3,YP,params);

    % theta_final = mod(theta_final,2*pi);
    % if theta_final<0
    %     theta_final = theta_final + 2*pi;
    % end

    % Valori target
    final_time_target = 2.9;
    final_time = sum(YP(1:3));
    r_target = 1.2;
    vt_target = sqrt(params.mu / r_target);
    vr_target = 0;
    theta_target = pi + sqrt(params.mu / r_target^3)*final_time;
    % disp(theta_target)
    % disp(theta_final)
    er_theta = mod(theta_final - theta_target + pi, 2*pi) - pi;
    % er_theta = mod(er_theta,2*pi);
    % if theta_target<0
    %     theta_target = theta_target + 2*pi;
    % end

    % Calcolo errori
    ER = [% H12;
          % H23;
          swf1;
          swf2;
          % H3final;
          final_time - final_time_target;
          r_final - r_target;
          er_theta;
          % lambda_theta_final;
          vr_final - vr_target;
          vt_final - vt_target;
          % lambda_m_final - 1;
          lambda_norm_sq - 1
          ];

    % %  Re-scaling
    % L = [  0.1;    % H12 ~ O(0.1)
    %        1;    % H23 ~ O(1)
    %        0.1;    % time error ~ O(0.1)
    %        1;    % radial pos error ~ O(1)
    %        0.1;    % radial vel error ~ O(0.1)
    %        0.01;    % tangential vel error ~ O(0.01)
    %        0.1;    % λm error ~ O(0.1)
    %        0.1];   % θ error ~ O(0.1 rad)
    % 
    % W = 1./L;
    % ER = W .* ER;
end

function s_right = discon(s_left, YP, icmp, params)
    % Nessuna discontinuità (stato continua tra gli archi)
    s_right = s_left;
end

function post_process(prb, sol, params)
% S = [r, theta, vr, vt, m, λr, λtheta, λvr, λvt, λm]
%     [1,     2,  3,  4, 5,  6,      7,   8,   9, 10]

%% Results analysis and visualization

% Extract solution data
ER = sol.ER;
YP = sol.YP;

% Calculate final mass
mass_3final = sol.S{end}(end,5);     % Mass at end


% %% Display results
% %YP=[tau1,tau2,tau3,λr(0),λvr(0),λvt(0),λm(0),λtheta(0)]
% 
% fprintf('\n--- Optimization Results ---\n');
% fprintf('Final parameters (YP):\n'); 
% YP'
% ER'
% fprintf('tau1: %12.4e\ntau2: %12.4e\ntau3: %12.4e\nλ_r(0): %12.4e\nλ_theta(0): %12.4e\nλ_vr(0): %12.4e\nλ_vt(0): %12.4e\nλ_m(0): %12.4e\n\n', YP);
% fprintf('Boundary condition residuals (ER):\n');
% fprintf('H12: %12.4e\nH23: %12.4e\nFinal Time: %12.4e\nPosition: %12.4e\nTheta costate: %12.4e\nRadial velocity: %12.4e\nTangential velocity: %12.4e\nMass costate: %12.4e\n\n', ER);
fprintf('\nFinal mass: %f\n', mass_3final);
%% Plotting
colors = {'b','r','g'};

%{ 
figure('Name','Trajectory Profile')
hold on; grid on; box on
for icmp = 1:prb.NCMP
    S = sol.S{icmp};
    plot(S(:,1), S(:,3), 'Color', colors{icmp}, 'LineWidth', 1.5)
end
xlabel('Radial Position (r)')
ylabel('Radial Velocity (v_r)')
title('Radial Phase Plane')
legend('1st Arc (Thrust)','2nd Arc (Coast)','3rd Arc (Thrust)')

figure('Name','Mass Evolution')
hold on; grid on; box on
for icmp = 1:prb.NCMP
    T = sol.T{icmp};
    S = sol.S{icmp};
    plot(T, S(:,5), 'Color', colors{icmp}, 'LineWidth', 1.5)
end
xlabel('Time')
ylabel('Mass (m)')
title('Mass Consumption History')


%}


figure('Name','Thrust Direction History')
hold on; grid on; box on
for icmp = 1:prb.NCMP
    T = sol.T{icmp};
    S = sol.S{icmp};
    phi = atan2(-S(:,9), -S(:,8));
    % phi = atan2(S(:,9), S(:,8));
    phi=unwrap(phi);
    plot(T, rad2deg(phi), 'Color', colors{icmp}, 'LineWidth', 1.5)
end
xlabel('Time')
ylabel('Thrust Angle (deg)')
title('Optimal Control History')



% Switching function
figure('Name','SWF')
hold on; grid on; box on

% % Plot trajectory arcs
% for icmp = 1:prb.NCMP
%     T = sol.T{icmp};
%     S = sol.S{icmp};
%     lambda_vr = S(:,8);  
%     lambda_vt = S(:,9);
%     % lambda_v = sqrt(lambda_vr.^2 + lambda_vt.^2);
%     m = S(:,5);
%     lambda_m = S(:,10);
%     C = params.C;
%     % swf = lambda_v./m - lambda_m./C;
%     swf = sqrt(lambda_vr.^2 + lambda_vt.^2)./m - lambda_m./C;
%     plot(T, swf, 'Color', colors{icmp}, 'LineWidth', 1.5)
% end

for icmp = 1:prb.NCMP
    T = sol.T{icmp};
    S = sol.S{icmp}; 
    swf = zeros(length(S),1);
    for i = 1:length(S)
        swf(i) = switchingF(1, S(i,:), icmp, YP, params);
    end
    plot(T, swf, 'Color', colors{icmp}, 'LineWidth', 1.5)
end
xlabel('T Time')
ylabel('swf')
title('Switching Function')
legend('1st Arc (Thrust)','2nd Arc (Coast)','3rd Arc (Thrust)')




% New 2D Trajectory Plot with Orbits
figure('Name','2D Trajectory with Orbits')
hold on; grid on; box on

% Plot trajectory arcs
for icmp = 1:prb.NCMP
    S = sol.S{icmp};
    theta = S(:,2);  % True anomaly from state
    r = S(:,1);      % Radial position
    [x,y] = pol2cart(theta, r);
    plot(x, y, 'Color', colors{icmp}, 'LineWidth', 1.5)
end

% Plot circular orbits
theta = linspace(0, 2*pi, 100);
% Departure orbit (r=1)
[x_dep, y_dep] = pol2cart(theta, ones(size(theta)));
plot(x_dep, y_dep, 'k--', 'LineWidth', 1)
% Arrival orbit (r=1.2)
[x_arr, y_arr] = pol2cart(theta, 1.2*ones(size(theta)));
plot(x_arr, y_arr, 'm--', 'LineWidth', 1)

xlabel('X Position')
ylabel('Y Position')
title('2D Trajectory with Circular Orbits')
legend('1st Arc (Thrust)','2nd Arc (Coast)','3rd Arc (Thrust)',...
       'Departure Orbit (r=1)', 'Arrival Orbit (r=1.2)')
axis equal
end