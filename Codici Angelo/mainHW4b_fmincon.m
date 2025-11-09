%% Transfer orbit optimization using fmincon (time-minimum problem)
clc;
close all;
clear;

%% Parameters

save_pic = false; 

DU = 6378.136;
TU = 806.81;

g0 = 9.8065e-3*(TU^2/DU);
par.c = 30*(TU/DU);
par.mu = 1;
par.a0    = 1e4/DU;

% g0 = 9.81;
% par.mu = 3.986e14;
% par.c = 30000;
% par.a0 = 10e6;
par.uTmax = 0.01 * g0;
par.e0    = 0.3;
par.lambda_theta = 0;
par.opz   = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Target values
% target.rf  = 42164e3;
target.rf = 42164/DU;
target.vrf = 0;
target.vtf = sqrt(par.mu/target.rf);
target.Hf  = -1;

% Initial guess for [theta0; lambda_vr0; lambda_vt0; tf]
YP_guess = [pi/2; 0; -1; 60];


%% fmincon setup
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'MaxIterations',300, ...
    'MaxFunctionEvaluations',1e9, ...
    'Algorithm','interior-point');

A = []; b = [];
Aeq = []; beq = [];
lb = []; ub = [];

% Run constrained optimization: minimize tf subject to terminal constraints
[YP, fval, exitflag] = fmincon(@(YP)objFun(YP), ...
                               YP_guess, ...
                               A, b, Aeq, beq, lb, ub, ...
                               @(YP)nonlcon(YP,par,target), ...
                               options);


%% Integrate the optimal solution

s0  = setIC(YP, par);
tf  = YP(4);
tspan = [0 tf];
[time, sol] = ode113(@(t,s) odeHW4(t, s, par), tspan, s0, par.opz);

% Extract states
r = sol(:,1);
thetastar = sol(:,2);
vr = sol(:,3);
vt = sol(:,4);
lambda_vr = sol(:,7);
lambda_vt = sol(:,8);
final_state = sol(end,:)';

fprintf('\n--- Optimized Parameters ---\n');
fprintf('Initial true anomaly θ*(0)          : %8.3f deg\n', rad2deg(YP(1)));
fprintf('Costate λ_vr(0) (radial velocity)   : %8.3e TU²/DU\n', YP(2));
fprintf('Costate λ_vt(0) (tang. velocity)    : %8.3e TU²/DU\n', YP(3));
fprintf('Optimized time of flight tf         : %8.3f Time Unit\n', YP(4));

%% Plotting 

% Creating folder for saving picutres
if ~exist('Figure','dir') && save_pic
    mkdir('Figure')
end

filename = fullfile('Figure\\fmincon_State.png');
fig=figure('name',filename,'NumberTitle', 'off');
subplot(2,2,1); plot(time, r); xlabel('time [TU]'); ylabel('r [DU]'); title('Radius');
subplot(2,2,2); plot(time, rad2deg(thetastar)); xlabel('time [TU]'); ylabel('\theta^* [deg]'); title('True Anomaly');
subplot(2,2,3); plot(time, vr); xlabel('time [TU]'); ylabel('v_r [DU/TU]'); title('Radial Velocity');
subplot(2,2,4); plot(time, vt); xlabel('time [TU]'); ylabel('v_t [DU/TU]'); title('Transverse Velocity');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end

filename = fullfile('Figure\\fmincon_StateErrors.png');
fig=figure('name',filename,'NumberTitle', 'off');
subplot(3,1,1); plot(time, r - target.rf); xlabel('time [TU]'); ylabel('r [DU]'); title('Radius - Target Radius');
subplot(3,1,2); plot(time, vr - target.vrf); xlabel('time [TU]'); ylabel('v_r [DU/TU]'); title('Radial Velocity - Target Radial Velocity');
subplot(3,1,3); plot(time, vt - target.vtf); xlabel('time [TU]'); ylabel('v_t [DU/TU]'); title('Transverse Velocity - Target Transverse Velocity');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end

filename = fullfile('Figure\\fmincon_Orbit.png');
fig=figure('name',filename,'NumberTitle', 'off');
[x, y] = pol2cart(thetastar, r);
plot(x, y, 'LineWidth', 1.5);
hold on;
theta_circle = linspace(0,2*pi,300);
p = par.a0*(1-par.e0^2);
r_init = p./(1+par.e0*cos(theta_circle));
[xi, yi] = pol2cart(theta_circle, r_init);
plot(xi, yi, 'k--', 'LineWidth', 1);
r_final = target.rf;
[xf, yf] = pol2cart(theta_circle, r_final*ones(size(theta_circle)));
plot(xf, yf, 'r--', 'LineWidth', 1);
axis equal;
xlabel('X Position [DU]'); ylabel('Y Position [DU]'); title('2D Transfer Trajectory');
legend('Transfer','Initial Orbit','Final Orbit');
hold off;
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end


%% Attitude: request (e)
phi = atan2(-lambda_vt, -lambda_vr);

phi=unwrap(phi);
quat = zeros(4, length(time));
bryant = zeros(3, length(time));

for i = 1:length(time)
    rotB = rot(3*pi/2, 1)*rot(phi(i), 3);
    rotNP = rot(thetastar(i), 3);
    rotNC = rotB * rotNP;
    quat(:,i) = dcm2quat(rotNC)';
    [a1, a2, a3] = dcm2angle(rotNC, 'ZYX'); 
    bryant(:,i) = [a1, a2, a3]'; 
    % bryant_lim=[pi;pi/2;pi];
    % for j=1:3
    %     while bryant(j,i) >= bryant_lim(j)
    %         bryant(j,i) = bryant(j,i) - 2*pi;
    %     end
    %     while bryant(j,i) < -bryant_lim(j)
    %         bryant(j,i) = bryant(j,i) + 2*pi;
    %     end
    % end
    % bryant(1,i)=unwrap(bryant(1,i));
end
for i=2:length(time)
    % Flip sign if discontinuity detected
    if dot(quat(:,i), quat(:,i-1)) < 0
        quat(:,i) = -quat(:,i);
    end
end
bryant=unwrap(bryant,[],2);
bryant = rad2deg(bryant);

%% Plotting: request (e)
filename = fullfile('Figure\\fmincon_quaternions.png');
fig=figure('name',filename,'NumberTitle', 'off');
subplot(2,2,1); plot(time, quat(1,:)); title('q_{0}'); xlabel('time [TU]'); ylabel('q_{0}');
subplot(2,2,2); plot(time, quat(2,:)); title('q_{1}'); xlabel('time [TU]'); ylabel('q_{1}');
subplot(2,2,3); plot(time, quat(3,:)); title('q_{2}'); xlabel('time [TU]'); ylabel('q_{2}');
subplot(2,2,4); plot(time, quat(4,:)); title('q_{3}'); xlabel('time [TU]'); ylabel('q_{3}');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end

filename = fullfile('Figure\\fmincon_Bryant.png');
fig=figure('name',filename,'NumberTitle', 'off');
plot(time, bryant); legend('Yaw (ψ)','Pitch (θ)','Roll (φ)'); title('Bryant Angles'); xlabel('time [TU]'); ylabel('angle [deg]');
% subplot(4,1,[1 2]); plot(time, bryant(1,:)); title('Yaw (ψ)'); xlabel('Time [s]'); ylabel('ψ [deg]');
% subplot(4,1,3); plot(time, bryant(2,:)); title('Pitch (θ)'); xlabel('Time [s]'); ylabel('θ [deg]');
% subplot(4,1,4); plot(time, bryant(3,:)); title('Roll (φ)'); xlabel('Time [s]'); ylabel('φ [deg]');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end

filename = fullfile('Figure\\fmincon_ThrustAngle.png');
fig=figure('name',filename,'NumberTitle', 'off');
plot(time,rad2deg(phi)); title('Thrust Angle'); xlabel('time [TU]'); ylabel('φ [deg]');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end


%% Function Definition

%objective function: final time to be minimized
function J = objFun(YP)
    J = YP(4);
end


% Nonlinear constraints: terminal boundary conditions
function [cin, ceq] = nonlcon(YP, par, target)
    s0 = setIC(YP, par);
    tf = YP(4);
    [~, sol] = ode113(@(t,s) odeHW4(t, s, par), [0 tf], s0, par.opz);
    final = sol(end,:)';
    lam=final(6:8);
    % Inequality constraints
    % cin = [hamiltonian(final, par)];
    cin=[];
    % Equality constraints: r, vr, vt, Hamiltonian
    ceq = [ final(1) - target.rf;           ... % radius
            final(3) - target.vrf;          ... % radial vel
            final(4) - target.vtf;          ... % transverse vel
            % hamiltonian(final, par) - target.Hf  % H(tf) = -1
            sum(lam.^2)-1;
            ];
end


% form the boundary condition from the given and the guessed boundary conditions
function s0 = setIC(YP, par)
    % initial state and costates
    a0 = par.a0; e0 = par.e0; mu = par.mu;
    p0 = a0 * (1 - e0^2);
    theta0 = YP(1);
    r0 = p0 / (1 + e0 * cos(theta0));
    vr0 = sqrt(mu/p0) * e0 * sin(theta0);
    vt0 = sqrt(mu/p0) * (1 + e0 * cos(theta0));
    xi0 = 1;
    
    lambda_vr0 = YP(2);
    lambda_vt0 = YP(3);
    lambda_r0 = ((1+e0*cos(theta0))^2 / (p0*sin(theta0)))*sqrt(mu/p0)*(lambda_vt0*sin(theta0) - lambda_vr0*cos(theta0));
    s0 = [r0; theta0; vr0; vt0; xi0; lambda_r0; lambda_vr0; lambda_vt0];
end


% ODEs
function dsdt = odeHW4(t, s, par)
    % State: [r; theta; vr; vt; xi]; 
    % Costate: [lambda_r; lambda_vr; lambda_vt]
    mu = par.mu;
    uTmax = par.uTmax;
    c = par.c;
    lambda_theta = par.lambda_theta;

    r = s(1); vr = s(3); vt = s(4); xi = s(5);
    lambda_r = s(6); lambda_vr = s(7); lambda_vt = s(8);
    

    % Control angle
    phi = atan2(-lambda_vt, -lambda_vr);

    % State dynamics
    drdt = vr;
    dthetadt = vt / r;
    dvrdt = vt^2 / r - mu / r^2 + uTmax / xi * cos(phi);
    dvtdt = -vr * vt / r + uTmax / xi * sin(phi);
    dxidt = - uTmax/c;

    % Costate dynamics
    dlambda_r = lambda_vr * (vt^2 - 2*mu/r) / r^2 - lambda_vt * (vr*vt) / r^2 + lambda_theta * vt / r^2;
    dlambda_vr = -lambda_r + lambda_vt * vt / r;
    dlambda_vt = -2 * lambda_vr * vt / r + lambda_vt * vr / r - lambda_theta / r;

    dsdt = [drdt; dthetadt; dvrdt; dvtdt; dxidt; dlambda_r; dlambda_vr; dlambda_vt];
end


% Hamilton's function
function H = hamiltonian(s, par)
    mu = par.mu; uTmax = par.uTmax; c=par.c;
    r = s(1); vr = s(3); vt = s(4); xi = s(5);
    lambda_r = s(6); lambda_vr = s(7); lambda_vt = s(8); 
    
    lam = [lambda_r; 0 ; lambda_vr; lambda_vt; 0];

    phi = atan2(-lambda_vt, -lambda_vr);
    f = [vr;
         vt/r;
         vt^2/r - mu/r^2 + uTmax/xi*cos(phi);
         -vr*vt/r + uTmax/xi*sin(phi);
         - uTmax/c];
        
    % lam = [lambda_r; 0 ; lambda_vr; lambda_vt];
    % 
    % phi = atan2(-lambda_vt, -lambda_vr);
    % f = [vr;
    %      vt/r;
    %      vt^2/r - mu/r^2 + uTmax/xi*cos(phi);
    %      -vr*vt/r + uTmax/xi*sin(phi)];

    H = lam' * f;
end


% Elementary Rotation Matrix
function rotM = rot(angle,axis)
    switch axis
        case 1
            rotM = [1 0 0;
                    0 cos(angle) sin(angle);
                    0 -sin(angle) cos(angle)];
        case 2
            rotM = [cos(angle) 0 -sin(angle);
                    0 1 0;
                    sin(angle) 0 cos(angle)];
        case 3
            rotM = [cos(angle) sin(angle) 0;
                    -sin(angle) cos(angle) 0;
                    0 0 1];
    end

end