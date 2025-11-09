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




function Energy_fuel_switch_system_RTN
clear
close all
clc

% save_pic = false;
params.PlotsOn = false;

% Possibly adding H3final into ER
prb.PBIS = 1;      % Bisection coefficient (>1)
prb.RMIN = 1e-2;    % Coefficient of Newton method: x_{k+1} = x_{x} - RMIN * f(x_k)/Df_{k}
prb.JMAX = 10;    % Maximum number of steps

prb.NY = 10;       % Number of differential equations
prb.KP = 1;        % Number of unknowns that are not initial boundary values (e.g., the time-length of an arc)
prb.K = 5 + 1;       % Number of unknowns (total). Must be equal to the number of boundary constraints
prb.NCMP = 1;      % Number of arcs (>=1)
prb.ATOL = 1e-10;   % Absolute tolerance for ODE solver
prb.RTOL = 1e-12;   % Relative tolerance for ODE solver

params.mu = 1;                    % Gravitational Parameters
params.Q = 0.1;                   % Max flux
params.C = 1;                     % Exhaust velocity
params.THR = params.Q*params.C;   % Thrust

prb.setIC = @setIC;
prb.funz = @funz;
prb.discon = @discon;
prb.bound = @bound;

prb.fsolve_opz = optimset('Display','Iter');
%fsolve_opz.Algorithm = 'levenberg-marquardt';
prb.fsolve_opz.MaxIter = 100;
prb.fsolve_opz.MaxFunEvals=1e6;
% prb.fsolve_opz.Display = 'final-detailed';
prb.solve = 1;

if exist('resultsRTN.mat','file')
    delete('resultsRTN.mat')
end
if exist('fuelOPT_resultsRTN.mat','file')
    delete('fuelOPT_resultsRTN.mat')
end

for t = 6:1:16
    params.TOF=t;
    TOF_i = (pi+sqrt(1/1.2^3)*params.TOF)/sqrt(params.mu/1^3);
    
    %[taui,...,λr(0),λtheta(0),λvr(0),vt(0),λxi(0)]
    prb.YP_guess = [TOF_i 0.01 0.01 0.01 0.01 0.01];
    
        
    for key = 0:0.05:1 % k=0 --> r_f=1 | k=1 --> r_f=1.2
        key
        params.epsl = 1;
        params.k=key;
        prb.params = params; 
        sol = solve_by_shooting(prb);
        prb.YP_guess=sol.YP;
        post_process(prb,sol, params)
    end
     
     
    
    for epslon = 1:-0.03:0.01 % eps=1 --> energy-optimal | eps=0 --> fuel-optimal
    % eps_vec = 1-linspace(0,1,15).^(2/3);
    % for epslon = eps_vec
        epslon
        params.epsl = epslon;
        params.k=1;
        prb.params = params;
        sol = solve_by_shooting(prb);
        prb.YP_guess=sol.YP;
        post_process(prb,sol, params)
    end
    
    for epslon = 0.01:-0.001:0.001 % eps=1 --> energy-optimal | eps=0 --> fuel-optimal
    % eps_vec = 1-linspace(0,1,15).^(2/3);
    % for epslon = eps_vec
        epslon
        params.epsl = epslon;
        params.k=1;
        prb.params = params;
        sol = solve_by_shooting(prb);
        prb.YP_guess=sol.YP;
        post_process(prb,sol, params)
    end
end



%% Problem dependant functions
% state = [r, th, vr, vt, xi, λr, λth, λvr, λvt, λxi]
%         [1, 2 , 3 , 4,  5,  6,  7,   8  , 9 ,  10]

    function [val, isExplicit] = setIC(j)
        % Set explicit value of some initial conidition s(j)
        % if s(j) == val, isExplicit must be true

        switch j
            case 1 % r
                isExplicit = 1; val = 1;
            case 2 % theta
                isExplicit = 1; val = 0;
            case 3 % vr
                isExplicit = 1; val = 0;
            case 4 % vt
                isExplicit = 1; val = 1;
            case 5 % xi
                isExplicit = 1; val = 1;
            otherwise
                val = NaN;
                isExplicit = false;
        end
    end

    function ds = funz(t, s, icmp, YP, params)
        % ODE function to be integrated.
        % state = [r, th, vr, vt, xi, λr, λth, λvr, λvt, λxi]
        %         [1, 2 , 3 , 4,  5,  6,  7,   8  , 9 ,  10]

        a_max = params.THR;
        C = params.C;
        mu = params.mu;
        epsl = params.epsl;

        tau = YP(icmp);

        % extract from state
        r = s(1);
        theta = s(2);
        vr = s(3);
        vt = s(4);
        xi = s(5);

        lambda_r = s(6);
        lambda_theta = s(7);
        lambda_vr = s(8);
        lambda_vt = s(9);
        lambda_xi = s(10);

        p = sqrt(lambda_vr^2+lambda_vt^2);
        ur = lambda_vr/p;
        ut = lambda_vt/p;
        swf = p/xi - lambda_xi/C - 1;

        if swf > epsl
            beta = 1;
        elseif swf < -epsl
            beta = 0;
        else
            beta = 0.5 * (1 + swf/epsl);
        end

        % state dynamics
        drdt = vr;
        dthetadt = vt/r;
        dvrdt = (vt^2)/r - mu/r^2 + a_max*beta*ur/xi;
        dvtdt = (-vr*vt)/r + a_max*beta*ut/xi;
        dxidt = -a_max*beta/C;

        % costates dynamics
        dlambda_r = -(2*lambda_vr*mu - lambda_vr*r*vt^2 - lambda_theta*r*vt + lambda_vt*r*vr*vt)/r^3;
        dlambda_theta = 0;
        dlambda_vr = (lambda_vt*vt)/r - lambda_r;
        dlambda_vt = -(lambda_theta - lambda_vt*vr + 2*lambda_vr*vt)/r;
        dlambda_xi = (a_max*beta*(lambda_vr*ur + lambda_vt*ut))/xi^2;
        

        ds = tau*[drdt; dthetadt; dvrdt; dvtdt; dxidt; dlambda_r; dlambda_theta; dlambda_vr; dlambda_vt; dlambda_xi];
    end



    function ER = bound(S_left, S_right, YP, params)
        % Computes the error between the final values (from integration)
        % and the final boundary condition.
        % state = [r, th, vr, vt, xi, λr, λth, λvr, λvt, λxi]
        %         [1, 2 , 3 , 4,  5,  6,  7,   8  , 9 ,  10]
        
        mu = params.mu;
        k = params.k;

        % Condizioni finali (t=7)
        s_final = S_right(:, end);
        r_final = s_final(1);
        theta_final = s_final(2);
        vr_final = s_final(3);
        vt_final = s_final(4);
        xi_final = s_final(5);
        lambda_r_final = s_final(6);
        lambda_theta_final = s_final(7);
        lambda_vr_final = s_final(8);
        lambda_vt_final = s_final(9);
        lambda_xi_final = s_final(10);
        % lambda_norm_sq = lambda_r_final^2+lambda_theta_final^2+lambda_vr_final^2+lambda_vt_final^2+lambda_m_final^2;
        % lambda0_norm_sq = sum(S_left(6:10,1).^2);
        tau=YP(1);

        TOF_target = params.TOF;

                            %[r theta vr vt     TOF    ]
        % final_state_target = [1.2 0 0 pi TOF_target];
        final_state_target = [1.2 pi+sqrt(mu/1.2^3)*TOF_target 0 sqrt(mu/1.2) TOF_target];

        TOF_initial = (pi+sqrt(mu/1.2^3)*TOF_target)/sqrt(mu/1^3);
        final_state_initial = [1 pi+sqrt(mu/1.2^3)*TOF_target 0 sqrt(mu/1) TOF_initial];

        iter_state = (1-k)*final_state_initial + k*final_state_target;

        ER = zeros([6,1]);
        ER(1) = tau - iter_state(5);
        ER(2) = r_final - iter_state(1);
        ER(3) = theta_final - iter_state(2);
        ER(4) = vr_final - iter_state(3);
        ER(5) = vt_final - iter_state(4);
        ER(6) = lambda_xi_final;
       
    end



    function s_right = discon(s_left, YP, icmp, params)
        s_right = s_left;
    end


    function post_process(prb, sol, params)
        % S = [r, theta, vr, vt, m, λr, λtheta, λvr, λvt, λm]
        %     [1,     2,  3,  4, 5,  6,      7,   8,   9, 10]

        %% Results analysis and visualization

        PlotsOn = params.PlotsOn;

        % Extract solution data
        ER = sol.ER;
        YP = sol.YP;
        mass_final = sol.S{end}(end,5);     % Mass at end

        YP
        fprintf('\nFinal mass: %f\n', mass_final);

        if exist('resultsRTN.mat','file')
            data = load('resultsRTN.mat');
            results = data.results;
        end

        TOF_name = sprintf('TOF%d',params.TOF);
        if (params.epsl - 0.001) < 1e-6
            if exist('fuelOPT_resultsRTN.mat','file')
                data = load('fuelOPT_resultsRTN.mat');
                fuel_results = data.fuel_results;
            end
            fuel_results.(TOF_name).final_mass=mass_final;
            fuel_results.(TOF_name).YP=YP;
            fuel_results.(TOF_name).TOF=params.TOF;
            save("fuelOPT_resultsRTN.mat","fuel_results")
        end
            


        itername = sprintf('eps%d_k%d', round(params.epsl*1000), round(params.k*1000));
        results.(TOF_name).(itername).epsl=params.epsl;
        results.(TOF_name).(itername).k = params.k;
        results.(TOF_name).(itername).TOF = params.TOF;
        results.(TOF_name).(itername).finalMass = mass_final;
        results.(TOF_name).(itername).YP=YP;
       

        %% Plotting

                %{
        if arc1_isThrust
            colors = {'r','k'};
            sizes = {1.5, 1.5};
        else
            colors = {'k','r'};
            sizes = {1.5, 1.5};
        end

        if ~exist('Figure','dir') && save_pic
            mkdir('Figure')
        end


figure('Name','Thrust Direction History')
hold on; grid on; box on
for icmp = 1:prb.NCMP
    T = sol.T{icmp};
    S = sol.S{icmp};
    % phi = atan2(-S(:,9), -S(:,8));
    phi = atan2(S(:,8), S(:,9));
    plot(T, rad2deg(phi))
end
xlabel('Time')
ylabel('Thrust Angle (deg)')
title('Optimal Control History')
        %}


        if PlotsOn
            figure; hold on; grid on; box on
            xlabel('t [s]'); ylabel('\\beta [adim]');
            title(sprintf('Andamento di \\beta per k=%.2f e \\epsilon = %.3f',params.k ,params.epsl));
        end

        epsl = params.epsl;
        C = params.C;

        for icmp = 1:prb.NCMP
            S = sol.S{icmp};
            T = (sol.T{icmp}-(icmp-1))*YP(icmp) + sum(YP(1:icmp-1));

            % Stati
            r  = S(:,1);
            xi = S(:,5);

            % Costati
            lambda_vr = S(:,8);
            lambda_vt  = S(:,9);
            lambda_xi = S(:,10);

            p = sqrt(lambda_vr.^2+lambda_vt.^2);
            swf = p./xi - lambda_xi/C - 1;
            beta_NOsat = 0.5 * (1 + swf./epsl);
            beta = zeros(size(p));
            for time = 1:length(r)
                if swf(time) > epsl
                    beta(time) = 1;
                elseif swf(time) < -epsl
                    beta(time) = 0;
                else
                    beta(time) =  beta_NOsat(time);
                end
            end

            if PlotsOn
                % Plottaggio
                plot(T, beta, 'LineWidth', 1.5);
                hold on
                plot(T, beta_NOsat,'--', 'LineWidth', 1);
            end
            results.(TOF_name).(itername).beta=beta;
            results.(TOF_name).(itername).time=T;
        end




        if PlotsOn
            figure; hold on; grid on; box on
            xlabel('x [km]'); ylabel('y [km]'); axis equal
            title(sprintf('Oribta per k=%.2f e \\epsilon = %.3f',params.k ,params.epsl));
        end

        % Plot trajectory arcs
        for icmp = 1:prb.NCMP
            S = sol.S{icmp};
            theta = S(:,2);  % True anomaly from state
            r = S(:,1);      % Radial position
            [x,y] = pol2cart(theta, r);
            results.(TOF_name).(itername).coor = [x,y];
            if PlotsOn
                plot(x, y, 'LineWidth', 1.5)
            end
        end

        if PlotsOn
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
            % legend('1st Arc (Thrust)','2nd Arc (Coast)','3rd Arc (Thrust)', 'Departure Orbit (r=1)', 'Arrival Orbit (r=1.2)')
            axis equal
        end

        save("resultsRTN.mat","results")


    end

end