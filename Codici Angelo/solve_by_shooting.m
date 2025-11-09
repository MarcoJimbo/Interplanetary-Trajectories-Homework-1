function sol = solve_by_shooting(prb)
% User-defined constants and functions
% NY = number of differential equations
% KP = number of unknown parameters (not initial conditions)
% K = total number of unknowns (parameters + initial conditions)
% 
% params = structure containing all problem data that are used in the
% USER_DEFINED function
% 
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
% else, set  isExplicit = 0
% 
% Returns a structure
% sol.T{icmp} vector of normalized Times
% sol.S{icmp} vector of normalized Times
% 
% function ds = funz(t, s, icmp, YP, params)
% %%% Inputs
% % t (1) independent variable \in [icmp-1, icmp]
% % s (NY,1) dependent variable (typically s=[x;lambda]
% % icmp (1) arc index
% % YP (K,1) vector of unknowns
% % params =  user-defined structure
% %%% Outputs
% % ds (NY,1) = rhs of the differential equations 
% 
% 
% function ER = bound(S_left, S_right, YP, params)
% Compute the violation of the boundary conditions
% %%% Inputs
% S_left (NYxNCMP) values of the state vector at the left boundary of each arc
% S_right (NYxNCMP) values of the state vector at the right bundary of each arc
% % YP (K,1) vector of unknowns
% % params =  user-defined structure
% %%% Outputs
% ER (K,1) =  errors on the boundary conditions
% 
%
% function s_right = discon(s_left, YP, icmp, params)
% % compute the jumps across each boundary.
% %%% Inputs
% % s_right (nx1): value of the state before discontinuity at the boundary
% % icmp-1 (before the start of the icmp-th arc)
% %%% Outputs
% % s_right (nx1): value of the state after discontinuity at the boundary 
% icmp-1 (immediately after the start of the icmp-th arc)
% 
% 

RMIN = prb.RMIN;
PBIS = prb.PBIS;
JMAX = prb.JMAX;

% global NY KP K NCMP

params = prb.params;
NY = prb.NY;
KP = prb.KP;
K = prb.K;
NCMP = prb.NCMP;
ATOL = prb.ATOL;
RTOL = prb.RTOL;
YP_guess = prb.YP_guess;
if isempty(YP_guess)
    error('YP_guess is not defined')
end


if prb.solve == 1
    % YP = fsolve(@(YP)shooting(YP,prb), YP_guess, prb.fsolve_opz); %solve the root-finding problem

    [YP, ierr] = newton_raphson(@(YP)shooting(YP,prb), YP_guess, RMIN, PBIS, JMAX);
    fprintf('ierr = %d', ierr)
else
    YP= YP_guess;
end



%% Create final solution
opz = odeset('RelTol',ATOL, 'AbsTol',RTOL);

% Define the initial condition
s0_minus = eval_IC(YP, prb);

S_left = zeros(NY, NCMP);
S_right = zeros(NY, NCMP);

for icmp = 1:NCMP
    tspan = [icmp-1, icmp];
    s0_plus = prb.discon(s0_minus,YP, icmp, params);
    [T_icmp, S_icmp] = ode113(@(t,s)prb.funz(t,s,icmp,YP,params), tspan, s0_plus, opz);   % S(i, j) = state s(j) at time T(i)

    %save
    S_left(:,icmp) = S_icmp(1,:)';
    S_right(:,icmp) = S_icmp(end,:)';    
    sol.T{icmp} = T_icmp;
    sol.S{icmp} = S_icmp;

    % update
    s0_minus = S_icmp(end,:)';
end
ER = prb.bound(S_left, S_right, YP, params);
sol.ER = ER;
sol.YP = YP;

end




function ER = shooting(YP, prb)
NCMP = prb.NCMP;
params = prb.params;

opz = odeset('RelTol',1e-8, 'AbsTol',1e-8);

% Define the initial condition
s0_minus = eval_IC(YP, prb);

S_left = zeros(prb.NY, prb.NCMP);
S_right = zeros(prb.NY, prb.NCMP);

for icmp = 1:NCMP
    tspan = [icmp-1, icmp];
    s0_plus = prb.discon(s0_minus,YP, icmp, params);
    [T_icmp, S_icmp] = ode113(@(t,s)prb.funz(t,s,icmp,YP, params), tspan, s0_plus, opz);   % S(i, j) = state s(j) at time T(i)

    %save
    S_left(:,icmp) = S_icmp(1,:)';
    S_right(:,icmp) = S_icmp(end,:)';

    % update
    s0_minus = S_icmp(end,:)';
end
ER = prb.bound(S_left, S_right, YP, params);
end

function s0 = eval_IC(YP, prb)
s0 = zeros(prb.NY,1);
kk = prb.KP;
for j=1:prb.NY
    
    [val, isExplicit] = prb.setIC(j);
    if isExplicit
        s0(j) = val;
    else
        kk = kk+1;
        s0(j) = YP(kk);
    end

end

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
    tol = 1e-5;
    f0 = funz(x0);
    norm_f0 = norm(f0);
    ierr = 0;  % Initialize with success code
    
    for iter = 1:JMAX
        fprintf('Iter: %d\n', iter);
        fprintf('YP: '); fprintf('  %f', x0); fprintf('\n');
        fprintf('ER: '); fprintf('  %f', f0); fprintf('\n');


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
                fprintf('Bisection: |f_new| = %f, f_new: ', norm(f_new));
                fprintf('%f ', f_new);
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