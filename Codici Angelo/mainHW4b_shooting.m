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

PBIS = 2;       % Bisection coefficient (>1) 
RMIN = 0.2;    % Coefficient of Newton method: x_{k+1} = x_{x} - RMIN * f(x_k)/Df_{k} 
JMAX = 300;     % Maximum number of steps 

%% Initial guess for [theta0; lambda_vr0; lambda_vt0; tf]
YP_guess = [pi/2; 0; -1; 60];

%% Main solve
fsolve_opz = optimset('Display','Iter'); 
fsolve_opz.MaxIter = 1000; 
fsolve_opz.MaxFunEvals=1e6;
fsolve_opz.Algorithm = 'levenberg-marquardt';
% YP = fsolve(@(YP)shooting(YP,par,target), YP_guess, fsolve_opz); %solve the root-finding problem
[YP, ierr] = newton_raphson(@(YP)shooting(YP,par,target), YP_guess, RMIN, PBIS, JMAX);
 
%% Integrate solution
s0 = setIC(YP, par);
tf = YP(4);
tspan = [0 tf];
[time, sol] = ode113(@(t,s) odeHW4(t, s, par), tspan, s0, par.opz);
r = sol(:,1);
thetastar = sol(:,2);
vr = sol(:,3);
vt = sol(:,4);
lambda_vr = sol(:,7);
lambda_vt = sol(:,8);
final=sol(end,:);

fprintf('\n--- Optimized Parameters ---\n');
fprintf('Initial true anomaly θ*(0)          : %8.3f deg\n', rad2deg(YP(1)));
fprintf('Costate λ_vr(0) (radial velocity)   : %8.3e TU²/DU\n', YP(2));
fprintf('Costate λ_vt(0) (tang. velocity)    : %8.3e TU²/DU\n', YP(3));
fprintf('Optimized time of flight tf         : %8.3f Time Unit\n', YP(4));

% Plot state histories
% Creating folder for saving picutres
if ~exist('Figure','dir') && save_pic
    mkdir('Figure')
end


filename = fullfile('Figure\\shooting_State.png');
fig=figure('name',filename,'NumberTitle', 'off');
subplot(2,2,1); plot(time, r); xlabel('time [TU]'); ylabel('r [DU]'); title('Radius');
subplot(2,2,2); plot(time, rad2deg(thetastar)); xlabel('time [TU]'); ylabel('\theta^* [deg]'); title('True Anomaly');
subplot(2,2,3); plot(time, vr); xlabel('time [TU]'); ylabel('v_r [DU/TU]'); title('Radial Velocity');
subplot(2,2,4); plot(time, vt); xlabel('time [TU]'); ylabel('v_t [DU/TU]'); title('Transverse Velocity');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end


filename = fullfile('Figure\\shooting_StateErrors.png');
fig=figure('name',filename,'NumberTitle', 'off');
subplot(3,1,1); plot(time, r - target.rf); xlabel('time [TU]'); ylabel('r [DU]'); title('Radius - Target Radius');
subplot(3,1,2); plot(time, vr - target.vrf); xlabel('time [TU]'); ylabel('v_r [DU/TU]'); title('Radial Velocity - Target Radial Velocity');
subplot(3,1,3); plot(time, vt - target.vtf); xlabel('time [TU]'); ylabel('v_t [DU/TU]'); title('Transverse Velocity - Target Transverse Velocity');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end


% Plot 2D trajectory with initial & final orbits
filename = fullfile('Figure\\shooting_Trajectory.png');
fig=figure('name',filename,'NumberTitle', 'off');
[x, y] = pol2cart(thetastar, r);
plot(x, y, 'LineWidth', 1.5);
hold on;
% initial elliptic orbit
theta_circle = linspace(0, 2*pi, 300);
p = par.a0 * (1 - par.e0^2);
r_init = p ./ (1 + par.e0 * cos(theta_circle));
[xi, yi] = pol2cart(theta_circle, r_init);
plot(xi, yi, 'k--', 'LineWidth', 1);
% final circular orbit
r_final = target.rf;
[xf, yf] = pol2cart(theta_circle, r_final * ones(size(theta_circle)));
plot(xf, yf, 'r--', 'LineWidth', 1);
xlabel('X Position [DU]');
ylabel('Y Position [DU]');
title('2D Trajectory with Initial and Final Orbits');
axis equal;
legend('Transfer Trajectory', 'Initial Orbit', 'Final Orbit');
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
end
for i=2:length(time)
    % Flip sign if discontinuity detected
    if dot(quat(:,i), quat(:,i-1)) < 0
        quat(:,i) = -quat(:,i);
    end
end
bryant = unwrap(bryant,[],2);
bryant = rad2deg(bryant);



%% Plotting: request (e)
filename = fullfile('Figure\\shooting_quaternions.png');
fig=figure('name',filename,'NumberTitle', 'off');
subplot(2,2,1); plot(time, quat(1,:)); title('q_{0}'); xlabel('time [TU]'); ylabel('q_{0}');
subplot(2,2,2); plot(time, quat(2,:)); title('q_{1}'); xlabel('time [TU]'); ylabel('q_{1}');
subplot(2,2,3); plot(time, quat(3,:)); title('q_{2}'); xlabel('time [TU]'); ylabel('q_{2}');
subplot(2,2,4); plot(time, quat(4,:)); title('q_{3}'); xlabel('time [TU]'); ylabel('q_{3}');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end

filename = fullfile('Figure\\shooting_Bryant.png');
fig=figure('name',filename,'NumberTitle', 'off');
plot(time, bryant); legend('Yaw (ψ)','Pitch (θ)','Roll (φ)'); title('Bryant Angles'); xlabel('time [TU]'); ylabel('angle [deg]');
% subplot(4,1,[1 2]); plot(time, bryant(1,:)); title('Yaw (ψ)'); xlabel('Time [s]'); ylabel('ψ [deg]');
% subplot(4,1,3); plot(time, bryant(2,:)); title('Pitch (θ)'); xlabel('Time [s]'); ylabel('θ [deg]');
% subplot(4,1,4); plot(time, bryant(3,:)); title('Roll (φ)'); xlabel('Time [s]'); ylabel('φ [deg]');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end

filename = fullfile('Figure\\shooting_ThrustAngle.png');
fig=figure('name',filename,'NumberTitle', 'off');
plot(time,rad2deg(phi)); title('Thrust Angle'); xlabel('time [TU]'); ylabel('φ [deg]');
if save_pic
    exportgraphics(fig, filename, 'Resolution', 500);
end

%% Functions

function ER = shooting(YP, par, target)
    % Unpack and integrate
    s0 = setIC(YP, par);
    tf = YP(4);
    tspan = [0 tf];
    [~, sol] = ode113(@(t,s) odeHW4(t, s, par), tspan, s0, par.opz);
    final = sol(end, :)';
    rf = final(1);
    vrf = final(3);
    vtf = final(4);
    lam=final(6:8);
    Hf = hamiltonian(final, par);

    ER = [
        rf  - target.rf;
        vrf - target.vrf;
        vtf - target.vtf;
        % Hf  - target.Hf
        sum(lam.^2) - 1;
    ];
end

function s0 = setIC(YP, par)
    % initial state and costates
    a0 = par.a0; e0 = par.e0; mu = par.mu;
    p0 = a0 * (1 - e0^2);
    theta0 = YP(1);
    r0 = p0 / (1 + e0 * cos(theta0));
    vr0 = sqrt(mu/p0) * e0 * sin(theta0);
    vt0 = sqrt(mu/p0) * (1 + e0 * cos(theta0));
    eps0 = 1;
    
    lambda_vr0 = YP(2);
    lambda_vt0 = YP(3);
    lambda_r0 = ((1+e0*cos(theta0))^2 / (p0*sin(theta0)))*sqrt(mu/p0)*(lambda_vt0*sin(theta0) - lambda_vr0*cos(theta0));

    s0 = [r0; theta0; vr0; vt0; eps0; lambda_r0; lambda_vr0; lambda_vt0];
end

function dsdt = odeHW4(t, s, par)
    % State: [r; theta; vr; vt]; Costate: [lambda_r; lambda_theta; lambda_vr; lambda_vt]
    mu = par.mu;
    uTmax = par.uTmax;
    c = par.c;

    r = s(1); theta = s(2); vr = s(3); vt = s(4); eps = s(5);
    lambda_r = s(6); lambda_theta=par.lambda_theta;
    lambda_vr = s(7); lambda_vt = s(8);
    

    % Control angle
    phi = atan2(-lambda_vt, -lambda_vr);

    % State dynamics
    drdt = vr;
    dthetadt = vt / r;
    dvrdt = vt^2 / r - mu / r^2 + uTmax / eps * cos(phi);
    dvtdt = -vr * vt / r + uTmax / eps * sin(phi);
    depsdt = - uTmax/c;

    % Costate dynamics
    dlambda_r = lambda_vr * (vt^2 - 2*mu/r) / r^2 - lambda_vt * (vr*vt) / r^2 + lambda_theta * vt / r^2;
    dlambda_vr = -lambda_r + lambda_vt * vt / r;
    dlambda_vt = -2 * lambda_vr * vt / r + lambda_vt * vr / r - lambda_theta / r;

    dsdt = [drdt; dthetadt; dvrdt; dvtdt; depsdt; dlambda_r; dlambda_vr; dlambda_vt];
end

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


function df = compute_gradient(funz, x0, f0)
    eps = 1e-5;
    n = length(x0);
    m = length(f0);  % Get the number of functions in funz
    df = zeros(m, n); % Jacobian matrix is m × n (functions × variables)
    for i = 1:n
        dx = zeros(n,1);
        dx(i) = eps;
        x = x0 + dx;
        f = funz(x);
        df(:,i) = (f - f0) / eps; % Store column-wise
    end

    % Check if Jacobian is ill-conditioned
    rcond_val = rcond(df);
    if rcond_val < eps
        warning('Jacobian is ill-conditioned. Reciprocal condition number: %e', rcond_val);
    end
end




function [x0, ierr] = newton_raphson(funz, x0, RMIN, PBIS, JMAX)
    tol = 1e-12;
    f0 = funz(x0);
    norm_f0 = norm(f0);
    ierr = 0;  % Initialize with success code
    
    for iter = 1:JMAX
        fprintf('Iter: %d\n', iter);
        fprintf('YP: '); fprintf('  %8.4e', x0); fprintf('\n');
        fprintf('ER: '); fprintf('  %8.4e', f0); fprintf('\n');


        % Compute Jacobian
        J = compute_gradient(funz, x0, f0);
        
        % Solve linear system for Newton step
        dx = -RMIN * (J\f0);
        
        x_new = x0 + dx;
        f_new = funz(x_new);
        norm_f_new = norm(f_new);
        
        % Backtracking line search
        for ibis = 1:10
            if norm_f_new >= PBIS * norm_f0
                % fprintf('Bisection: |f_new| = %8.4e, f_new: ', norm(f_new));
                % fprintf('%8.4e ', f_new);
                fprintf('\n');
                
                dx = dx/2;
                x_new = x0 + dx;
                f_new = funz(x_new);
                norm_f_new = norm(f_new);
            else
                x0 = x_new;
                f0 = f_new;
                norm_f0 = norm_f_new;
                break;
            end
        end
        
        % Check for backtracking failure
        if ibis > 10
            ierr = -3;
            break;
        end
        
        % Check for convergence
        if norm_f0 < tol
            ierr = 0;  % Success code
            break;
        end
    end
    
    % If max iterations reached without convergence
    if iter == JMAX && norm_f0 >= tol
        ierr = -1;
    end
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