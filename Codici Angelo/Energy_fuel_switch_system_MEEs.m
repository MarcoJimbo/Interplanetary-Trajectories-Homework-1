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




function Energy_fuel_switch_system_MEEs
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

if exist('resultsMEE.mat','file')
    delete('resultsMEE.mat')
end
if exist('fuelOPT_resultsMEE.mat','file')
    delete('fuelOPT_resultsMEE.mat')
end

for t = 5:1:17
    params.TOF=t;
    TOF_i = (pi+sqrt(1/1.2^3)*params.TOF)/sqrt(params.mu/1^3);

    %[taui,...,λp(0),λex(0),λey(0),λL(0),λxi(0)]
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

    % prb.YP_guess=[5.4000 -0.1770 -0.0466 0.1115 0.0211 -0.1535];


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
        params.epsl = epslon;
        params.k=1;
        prb.params = params;
        sol = solve_by_shooting(prb);
        prb.YP_guess=sol.YP;
        post_process(prb,sol, params)
    end
end


%% Problem dependant functions
% state = [p, ex, ey, L, xi, λp, λex, λey, λL, λxi]
%         [1, 2 , 3 , 4,  5,  6,  7,   8 , 9 ,  10]

    function [val, isExplicit] = setIC(j)
        % Set explicit value of some initial conidition s(j)
        % if s(j) == val, isExplicit must be true

        switch j
            case 1 % p,
                isExplicit = 1; val = 1;
            case 2 % ex
                isExplicit = 1; val = 0;
            case 3 % ey,
                isExplicit = 1; val = 0;
            case 4 % L
                isExplicit = 1; val = 0;
            case 5 % xi
                isExplicit = 1; val = 1;
            otherwise
                val = NaN;
                isExplicit = false;
        end
    end

    function ds = funz(t, s, icmp, YP, params)
        % ODE function to be integrated.
        % state = [p, ex, ey, L, xi, λp, λex, λey, λL, λxi]
        %         [1, 2 , 3 , 4,  5,  6,  7,   8 , 9 ,  10]

        THR = params.THR;
        C = params.C;
        mu = params.mu;
        epsl = params.epsl;

        tau = YP(icmp);

        % extract from state
        p = s(1);
        ex = s(2);
        ey = s(3);
        L = s(4);
        xi = s(5);

        lambda_p = s(6);
        lambda_ex = s(7);
        lambda_ey = s(8);
        lambda_L = s(9);
        lambda_xi = s(10);    

        W = 1 + ex*cos(L) + ey * sin(L);
        nx = ex + cos(L);
        ny = ey + sin(L);
        lambda = [lambda_p lambda_ex lambda_ey lambda_L]';


        A = sqrt(mu/p)*[0 0 0 W^2/p]';

        B = sqrt(p/mu)*[   0        2*p/W  ;
                         sin(L) cos(L)+nx/W;
                        -cos(L) sin(L)+ny/W;
                           0          0    ];
        swf = - 1 + norm(B'*lambda)/xi - lambda_xi/C;
        a_max = THR; %DUBBIOOOOOOOOOOOO

        if swf > epsl
            beta = 1;
        elseif swf < -epsl
            beta = 0;
        else
            beta = 0.5 * (1 + swf/epsl);
        end

        u_vers = B'*lambda / norm(B'*lambda);

        ds = zeros([10,1]);
        % state dynamics
        ds(1:4) = A + B*u_vers*a_max*beta/xi; 
        ds(5) = -a_max*beta/C;
        
        % costates dynamics
        u1=u_vers(1);
        u2=u_vers(2);
        ds(6) = (lambda_L*(mu/p)^(1/2)*(ex*cos(L) + ey*sin(L) + 1)^2)/p^2 - (beta*lambda_ex*a_max*((u2*(cos(L) + (ex + cos(L))/(ex*cos(L) + ey*sin(L) + 1)))/(2*mu*(p/mu)^(1/2)) + (u1*sin(L))/(2*mu*(p/mu)^(1/2))))/xi + (beta*lambda_ey*a_max*((u1*cos(L))/(2*mu*(p/mu)^(1/2)) - (u2*(sin(L) + (ey + sin(L))/(ex*cos(L) + ey*sin(L) + 1)))/(2*mu*(p/mu)^(1/2))))/xi + (lambda_L*mu*(ex*cos(L) + ey*sin(L) + 1)^2)/(2*p^3*(mu/p)^(1/2)) - (2*beta*lambda_p*u2*a_max*(p/mu)^(1/2))/(xi*(ex*cos(L) + ey*sin(L) + 1)) - (beta*lambda_p*p*u2*a_max)/(mu*xi*(p/mu)^(1/2)*(ex*cos(L) + ey*sin(L) + 1));
        ds(7) = (2*beta*lambda_p*p*u2*a_max*cos(L)*(p/mu)^(1/2))/(xi*(ex*cos(L) + ey*sin(L) + 1)^2) - (2*lambda_L*cos(L)*(mu/p)^(1/2)*(ex*cos(L) + ey*sin(L) + 1))/p + (beta*lambda_ey*u2*a_max*cos(L)*(ey + sin(L))*(p/mu)^(1/2))/(xi*(ex*cos(L) + ey*sin(L) + 1)^2) - (beta*lambda_ex*u2*a_max*sin(L)*(ey + sin(L))*(p/mu)^(1/2))/(xi*(ex*cos(L) + ey*sin(L) + 1)^2);
        ds(8) = (2*beta*lambda_p*p*u2*a_max*sin(L)*(p/mu)^(1/2))/(xi*(ex*cos(L) + ey*sin(L) + 1)^2) - (beta*lambda_ey*u2*a_max*cos(L)*(ex + cos(L))*(p/mu)^(1/2))/(xi*(ex*cos(L) + ey*sin(L) + 1)^2) - (2*lambda_L*sin(L)*(mu/p)^(1/2)*(ex*cos(L) + ey*sin(L) + 1))/p + (beta*lambda_ex*u2*a_max*sin(L)*(ex + cos(L))*(p/mu)^(1/2))/(xi*(ex*cos(L) + ey*sin(L) + 1)^2);
        ds(9) = (2*beta*lambda_p*p*u2*a_max*(ey*cos(L) - ex*sin(L))*(p/mu)^(1/2))/(xi*(ex*cos(L) + ey*sin(L) + 1)^2) - (beta*lambda_ey*a_max*(u1*sin(L)*(p/mu)^(1/2) + u2*(p/mu)^(1/2)*(cos(L) + cos(L)/(ex*cos(L) + ey*sin(L) + 1) - ((ey*cos(L) - ex*sin(L))*(ey + sin(L)))/(ex*cos(L) + ey*sin(L) + 1)^2)))/xi - (2*lambda_L*(ey*cos(L) - ex*sin(L))*(mu/p)^(1/2)*(ex*cos(L) + ey*sin(L) + 1))/p - (beta*lambda_ex*a_max*(u1*cos(L)*(p/mu)^(1/2) - u2*(p/mu)^(1/2)*(sin(L) + sin(L)/(ex*cos(L) + ey*sin(L) + 1) + ((ey*cos(L) - ex*sin(L))*(ex + cos(L)))/(ex*cos(L) + ey*sin(L) + 1)^2)))/xi;
        ds(10) = (beta*a_max*(p/mu)^(1/2)*(lambda_ex*u1*sin(L) + 2*lambda_ey*u2*sin(L) + ex*lambda_ex*u2 + ey*lambda_ey*u2 + 2*lambda_p*p*u2 + 2*lambda_ex*u2*cos(L) - lambda_ey*u1*cos(L) + ex*lambda_ex*u2*cos(L)^2 - ex*lambda_ey*u1*cos(L)^2 + (ex*lambda_ex*u1*sin(2*L))/2 + (ex*lambda_ey*u2*sin(2*L))/2 + (ey*lambda_ex*u2*sin(2*L))/2 - (ey*lambda_ey*u1*sin(2*L))/2 + ey*lambda_ex*u1*sin(L)^2 + ey*lambda_ey*u2*sin(L)^2))/(xi^2*(ex*cos(L) + ey*sin(L) + 1));

        ds = tau*ds;
    end



    function ER = bound(S_left, S_right, YP, params)
        % Computes the error between the final values (from integration)
        % and the final boundary condition.
        % state = [p, ex, ey, L, xi, λp, λex, λey, λL, λxi]
        %         [1, 2 , 3 , 4,  5,  6,  7,   8 , 9 ,  10]
        
        mu = params.mu;
        k = params.k;

        % Condizioni finali (t=7)
        s_final = S_right(:, end);
        p_final = s_final(1);
        ex_final = s_final(2);
        ey_final = s_final(3);
        L_final = s_final(4);
        xi_final = s_final(5);
        lambda_p_final = s_final(6);
        lambda_ex_final = s_final(7);
        lambda_ey_final = s_final(8);
        lambda_L_final = s_final(9);
        lambda_xi_final = s_final(10);
        % lambda_norm_sq = lambda_r_final^2+lambda_theta_final^2+lambda_vr_final^2+lambda_vt_final^2+lambda_m_final^2;
        % lambda0_norm_sq = sum(S_left(6:10,1).^2);
        tau=YP(1);

        TOF_target = params.TOF;

                            %[p ex ey L     TOF    ]
        % final_state_target = [1.2 0 0 pi TOF_target];
        final_state_target = [1.2 0 0 pi+sqrt(mu/1.2^3)*TOF_target TOF_target];

        TOF_initial = (pi+sqrt(1/1.2^3)*TOF_target)/sqrt(mu/1^3);
        final_state_initial = [1 0 0 pi+sqrt(mu/1.2^3)*TOF_target TOF_initial];

        iter_state = (1-k)*final_state_initial + k*final_state_target;

        ER = zeros([6,1]);
        ER(1) = tau - iter_state(5);
        ER(2) = p_final - iter_state(1);
        ER(3) = ex_final - iter_state(2);
        ER(4) = ey_final - iter_state(3);
        ER(5) = L_final - iter_state(4);
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

        if exist('resultsMEE.mat','file')
            data = load('resultsMEE.mat');
            results = data.results;
        end

      TOF_name = sprintf('TOF%d',params.TOF);
        if (params.epsl - 0.001) < 1e-6
            if exist('fuelOPT_resultsMEE.mat','file')
                data = load('fuelOPT_resultsMEE.mat');
                fuel_results = data.fuel_results;
            end
            fuel_results.(TOF_name).final_mass=mass_final;
            fuel_results.(TOF_name).YP=YP;
            fuel_results.(TOF_name).TOF=params.TOF;
            save("fuelOPT_resultsMEE.mat","fuel_results")
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

        mu   = params.mu;
        THR  = params.THR;
        epsl = params.epsl;
        C    = params.C;

        for icmp = 1:prb.NCMP
            S = sol.S{icmp};
            T = (sol.T{icmp}-(icmp-1))*YP(icmp) + sum(YP(1:icmp-1));

            % Stati
            p  = S(:,1);
            ex = S(:,2);
            ey = S(:,3);
            L  = S(:,4);
            xi = S(:,5);

            % Costati
            lambda_p  = S(:,6);
            lambda_ex = S(:,7);
            lambda_ey = S(:,8);
            lambda_L  = S(:,9);
            lambda_xi = S(:,10);

            % Calcolo W, nx, ny
            W  = 1 + ex.*cos(L) + ey.*sin(L);
            nx = ex + cos(L);
            ny = ey + sin(L);

            % Matrici B e vettore lambda (righe → colonne)
            lambda = [lambda_p, lambda_ex, lambda_ey, lambda_L]';

            % Loop sui punti (perché B cambia punto per punto)
            beta = zeros(size(p));
            beta_NOsat = zeros(size(p));
            for time = 1:length(p)
                B = sqrt(p(time)/mu) * [  0         2*p(time)/W(time) ;
                    sin(L(time)) cos(L(time))+nx(time)/W(time);
                    -cos(L(time)) sin(L(time))+ny(time)/W(time);
                    0         0          ];
                lam = lambda(:,time);

                swf = - 1 + norm(B'*lam)/xi(time) - lambda_xi(time)/C;

                beta_NOsat(time) = 1/2 * (1+swf/epsl);

                % saturazione
                if swf > epsl
                    beta(time)=1;
                elseif swf<-epsl
                    beta(time) = 0;
                else
                    beta(time)=beta_NOsat(time);
                end
            end

            if PlotsOn
                % Plottaggio
                plot(T, beta, 'LineWidth', 1.5);
                hold on
                plot(T, beta_NOsat,'--', 'LineWidth', 1);
            end
            results.(itername).beta=beta;
            results.(itername).time=T;
        end


        if PlotsOn
            figure; hold on; grid on; box on
            xlabel('x [km]'); ylabel('y [km]'); axis equal
            title(sprintf('Oribta per k=%.2f e \\epsilon = %.3f',params.k ,params.epsl));
        end

        for icmp = 1:prb.NCMP
            S = sol.S{icmp};
            p  = S(:,1);
            ex = S(:,2);
            ey = S(:,3);
            L  = S(:,4);

            % eccentricità ed anomalia vera
            e     = sqrt(ex.^2 + ey.^2);
            a     = p./(1-e.^2);              %#ok<NASGU> % utile se vuoi
            argp  = 2*atan2(ey, ex+e);        % atan2 più robusto
            theta = L - argp;

            % distanza radiale
            r = p ./ (1 + e.*cos(theta));

            % coordinate cartesiane nel piano orbitale
            x = r .* cos(L);
            y = r .* sin(L);
            results.(itername).coor = [x,y];
            if PlotsOn
                % plot orbita
                plot(x, y, 'LineWidth', 1.5);
                hold on
                % Plot circular orbits
                theta = linspace(0, 2*pi, 100);
                % Departure orbit (r=1)
                [x_dep, y_dep] = pol2cart(theta, ones(size(theta)));
                plot(x_dep, y_dep, 'k--', 'LineWidth', 1)
                % Arrival orbit (r=1.2)
                [x_arr, y_arr] = pol2cart(theta, 1.2*ones(size(theta)));
                plot(x_arr, y_arr, 'm--', 'LineWidth', 1)
            end
        end

        save("resultsMEE.mat","results")


    end

end