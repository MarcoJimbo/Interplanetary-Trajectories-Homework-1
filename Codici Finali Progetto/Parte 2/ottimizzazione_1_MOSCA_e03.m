%% OTTIMIZZAZIONE TRAIETTORIA SATURNO - ENCELADO
clc; clear variables; close all;

%% 1. DATI GLOBALI E FISICI
global gamma_max W R_start r_umbriel mu_ur g0



% --- Parametri Fisici ---
mu_ur = 5793939;                 % [km^3/s^2] Parametro gravitazionale Urano
r_ur = 25362;                     % [km] Raggio Urano
r_umbriel = 266000;             % [km] Raggio orbita Umbriel
g0 = 9.81;                       % [m/s^2] Accelerazione grav. terrestre al livello del mare

%% PARAMETRI MOTORE (Definizione F/m)
F_max = 2;          % [N] Spinta massima del motore
m0_start = 850;      % [kg] Massa iniziale STIMATA per l'ottimizzazione
Isp = 3000; % [s] Impulso specifico
W = g0*Isp*10^(-3); % [km/s] Velocità di scarico del motore

% Calcolo di gamma_max: (F/m0) convertito in km/s^2
% La conversione N -> (kg * km / s^2) si ottiene dividendo F per 1000
gamma_max_ms2 = F_max / m0_start;          % Risultato in m/s^2 (ca. 6.15e-4)
gamma_max = gamma_max_ms2 / 1000;          % Risultato corretto in km/s^2


% --- Definizione Orbita di Partenza (ESTERNA) ---
% Partiamo da un'orbita eccentrica più alta di quella di Umbriel
R_start = r_umbriel + 400000;   % Partenza a rp~616000 km

%% 2. CONDIZIONI AL CONTORNO (SCALING)
% Stato Iniziale (Start)
r0 = R_start;

e0=0.8; % eccentricità orbita start
a0=r0/(1-e0);

theta0 = 0; 
vr0 = 0;
vtheta0 = sqrt((mu_ur*(1+e0))/r0); 

% Stato Finale Target (Umbriel) - Usato solo per riferimento qui
rf_target = r_umbriel;

% normalizzazione
L  = R_start;
V  = sqrt(mu_ur*(1+0.8) / L);
T  = L / V;
A  = V / T;  

%% 3. IMPOSTAZIONE OTTIMIZZAZIONE (FIRST GUESS & BOUNDS)



% --- BOUNDS (LIMITI) AGGIORNATI ---
% LB = [P1, P3, P4, P5, t_fin]
LB = [-1, -0.5, -0.5, -2, 1e6];    
UB = [1, 0.5, 0.5, 2, 5e8];

% --- OPZIONI FMINCON ---
fminconopts = optimoptions("fmincon", ...
  "Algorithm","interior-point", ...
  "Display","iter", ...
  "FiniteDifferenceType","central", ...
  "StepTolerance",1e-16, ...
  "OptimalityTolerance",1e-10, ...
  "MaxFunctionEvaluations",50000, ...
  "MaxIterations",2000);


% PSO
opts = optimoptions('particleswarm', ...
  'SwarmSize', 1000, ...
  'MaxIterations', 100, ...
  'MaxStallIterations', 80, ...
  'FunctionTolerance', 1e-6, ...
  'HybridFcn',{@fmincon,fminconopts} ,...
  'Display','iter');
bestfval = inf;

[x,fval] = particleswarm(@(X) fun(X),5, LB, UB,opts);


[ ~ , F] = Vincoli(x);      % funzione che ritorna il vettore dei residui
disp(F)
disp(norm(F))


sol = x;



%% 4. PROPAGAZIONE SOLUZIONE OTTIMA
% Estrazione variabili ottimizzate
P1_0 = sol(1); P3_0 = sol(2); P4_0 = sol(3); P5_0 = sol(4);
t_fin = sol(5)/T; 

% Setup integrazione finale ad alta precisione
X0 = [1, 0, 0, 1, 1]';
P0 = [P1_0; 0; P3_0; P4_0; P5_0]; 
IN = [X0; P0];

options = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
[t, integrator_output] = ode113(@(t, state) Dyn_Eqs(t,state), [0, t_fin], IN, options);

% --- LOGICA DI STOP  ---
t = t * T;
day_stop=floor(t(end)/60/60/24);
hour_stop=floor((t(end)/60/60)-day_stop*24);
T_fin_fin=t(end);

fprintf('Target Raggiunto a t =%d day ,%d hour (su un totale stimato di %.2f s).\n',day_stop,hour_stop,t_fin*T);
fprintf('Troncamento della traiettoria in eccesso.\n');

% Taglio tutti i vettori a questo indice
t = t(1:end);
integrator_output = integrator_output(1:end, :);

% -------------------------------------------------

% Estrazione Dati per Plotting (Ora i vettori hanno la lunghezza giusta)
x1 = integrator_output(:,1)*L;
x2 = integrator_output(:,2);
x3 = integrator_output(:,3)*V;
x4 = integrator_output(:,4)*V;
x5 = integrator_output(:,5);
P1 = integrator_output(:,6)/L;
P3 = integrator_output(:,8)/V;
P4 = integrator_output(:,9)/V;
P5 = integrator_output(:,10);

%% 5. POST-PROCESSING (Calcolo Controllo e Angoli)
H_star = [];
CONTROL = [];
DELTA = [];     % Angolo di spinta locale
Switch_Function = [];

for i = 1:length(t)
    SF = sqrt(P3(i)^2 + P4(i)^2) / x5(i) + P5(i) / W ;
    Switch_Function = [Switch_Function; SF];
    
    if SF > 0
        % MOTORE ACCESO
        gamma = gamma_max;
        CONTROL = [CONTROL, gamma];
        % Calcolo angolo (atan2 per il quadrante corretto)
        su2 = -P3(i) / sqrt(P3(i)^2+P4(i)^2);
        cu2 = -P4(i) / sqrt(P3(i)^2+P4(i)^2);
        u2 = atan2(su2,cu2);
        DELTA = [DELTA, u2];
    else 
        % MOTORE SPENTO
        gamma = 0;
        CONTROL = [CONTROL, gamma];
        DELTA = [DELTA, NaN]; % NaN per non plottare
    end
    
    % Calcolo Hamiltoniana (Verifica numerica)
    H = P1(i)*x3(i) + P3(i)*(x4(i)^2/x1(i) - mu_ur/x1(i)^2) + ...
          P4(i)*(-x4(i)*x3(i)/x1(i)) - gamma*SF;
    H_star = [H_star; H]; 
end

%% 6. GRAFICI DEI RISULTATI

% --- PLOT 1: HAMILTONIANA ---
figure(1);
yline(0,'k-','LineWidth',2)
plot(t, H_star, '-', 'LineWidth', 2,'Color',[0.7 0.1 0]);hold on;
title('Hamiltonian (Constant)');
xlabel('Time [s]'); ylabel('H');
xlim([0,max(t)])
grid on;

%% --- PLOT 2: P5 ---
figure(2);
yline(0,'k-','LineWidth',2)
plot(t, P5, '-', 'LineWidth', 2,'Color',[0.7 0.1 0]);hold on;
title('Lagrange multiplicator of x_5 = m / m_0');
xlabel('Time [s]'); ylabel('P_5');
xlim([0,max(t)])
grid on;

%%
% --- PLOT 2: TRAIETTORIA (Advanced) ---
figure(3); hold on; grid on; axis equal;
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k'); % Sfondo nero

x_plot = x1 .* cos(x2);
y_plot = x1 .* sin(x2);
alpha_ang = linspace(0, 2*pi, 720);
for i=1:length(alpha_ang)
r_ellipse=a0*(1-e0^2)/(1+e0*cos(alpha_ang(i)));
xell(i)=r_ellipse*cos(alpha_ang(i));
yell(i)=r_ellipse*sin(alpha_ang(i));
end
% Disegno Corpi Celesti e Orbite
plot(0, 0, 'bo', 'MarkerFaceColor', [0.1 0.6 0.9], 'MarkerSize', 12, 'DisplayName', 'Urano');
plot(r_umbriel*cos(alpha_ang), r_umbriel*sin(alpha_ang), 'Color', [1 0.6 0.2], 'LineWidth', 2, 'DisplayName', 'Target: Umbriel');
plot(xell, yell, 'Color', [0.7 0.5 0.5], 'LineStyle', '--', 'DisplayName', 'Start: Esterno','LineWidth', 2);

% Logica Disegno Archi (Spinta vs Coasting)
thrust_idx = find(CONTROL > 0);
arrow_len = 100000; % Lunghezza frecce (regolabile)

% Disegno base Coasting (tutta la linea in ciano tratteggiato)
plot(x_plot, y_plot, '--', 'LineWidth', 1, 'DisplayName', 'Coasting','Color',[0.1 0.7 0.1]);

if ~isempty(thrust_idx)
    % Identifica segmenti di spinta
    diff_idx = diff(thrust_idx);
    break_points = find(diff_idx > 1);
    start_arcs = [thrust_idx(1), thrust_idx(break_points + 1)];
    end_arcs = [thrust_idx(break_points), thrust_idx(end)];
    
    for k = 1:length(start_arcs)
        idx_s = start_arcs(k);
        idx_e = end_arcs(k);
        idx_m = round((idx_s + idx_e)/2); % Punto medio
        
        % Sovrascrivi parte spinta con linea rossa
        plot(x_plot(idx_s:idx_e), y_plot(idx_s:idx_e), 'LineWidth', 2, 'Color',[0.4 0.4 0.4],'HandleVisibility', 'on');
        
        % Markers On/Off
        plot(x_plot(idx_s), y_plot(idx_s), '^','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor', [0 0.7 0], 'MarkerSize', 6, 'HandleVisibility', 'on');
        plot(x_plot(idx_e), y_plot(idx_e), 's','MarkerEdgeColor',[0.7 0 0], 'MarkerFaceColor',[0.7 0 0], 'MarkerSize', 6, 'HandleVisibility', 'on');
        
        % --- Frecce direzione spinta: ~10 per arco ---
        N_arrows = 10;
        idx_list = unique(round(linspace(idx_s, idx_e, min(N_arrows, idx_e-idx_s+1))));
        
        for ia = idx_list
            u2 = DELTA(ia);                 % u2: pitch, positivo clockwise, in [-pi,pi]
            if ~isnan(u2)
        
                % wrap coerente in [-pi,pi] anche senza toolbox
                u2 = atan2(sin(u2), cos(u2));
        
                % Angolo inerziale (CCW da +x):
                % tangenziale = theta + pi/2, e u2 positivo clockwise => sottraggo u2
                ang = x2(ia) + pi/2 - u2;
        
                dx = cos(ang) * arrow_len;
                dy = sin(ang) * arrow_len;
        
                % Per non riempire la legenda con 100 frecce:
                hv = 'off';
                if (k==1) && (ia==idx_list(1)), hv = 'on'; end
        
                quiver(x_plot(ia), y_plot(ia), dx, dy, 0, ...
                    'Color',[0.7 0 0.7], 'LineWidth', 1, 'MaxHeadSize', 0.5, ...
                    'HandleVisibility', hv);
            end
        end

    end
end
plot(x_plot(idx_e), y_plot(idx_e),'o','MarkerEdgeColor', [0 0.8 0.8], 'MarkerFaceColor',[0 0.9 0.9],'HandleVisibility', 'off');      
title('Optimized Trajectory', 'Color', 'k');
xlabel('x [km]', 'Color', 'k'); ylabel('y [km]', 'Color', 'k');
legend('Uranus', 'Umbriel Orbit', 'Start', 'Coasting','Thrust','On','Off','Direction Arrow', 'Location', 'bestoutside', 'TextColor', 'k');


% plot zoom su traiettoria per leggibilità
figure(4); hold on; grid on; axis equal;
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k'); % Sfondo nero

x_plot = x1 .* cos(x2);
y_plot = x1 .* sin(x2);
alpha_ang = linspace(0, 2*pi, 720);
for i=1:length(alpha_ang)
r_ellipse=a0*(1-e0^2)/(1+e0*cos(alpha_ang(i)));
xell(i)=r_ellipse*cos(alpha_ang(i));
yell(i)=r_ellipse*sin(alpha_ang(i));
end
% Disegno Corpi Celesti e Orbite
plot(0, 0, 'bo', 'MarkerFaceColor', [0.1 0.6 0.9], 'MarkerSize', 12, 'DisplayName', 'Urano');
plot(r_umbriel*cos(alpha_ang), r_umbriel*sin(alpha_ang), 'Color', [1 0.6 0.2], 'LineWidth', 2, 'DisplayName', 'Target: Umbriel');
plot(xell, yell, 'Color', [0.7 0.5 0.5], 'LineStyle', '--', 'DisplayName', 'Start: Esterno','LineWidth', 2);

% Logica Disegno Archi (Spinta vs Coasting)
thrust_idx = find(CONTROL > 0);
arrow_len = 100000; % Lunghezza frecce (regolabile)

% Disegno base Coasting (tutta la linea in ciano tratteggiato)
plot(x_plot, y_plot, '--', 'LineWidth', 1, 'DisplayName', 'Coasting','Color',[0.1 0.7 0.1]);

if ~isempty(thrust_idx)
    % Identifica segmenti di spinta
    diff_idx = diff(thrust_idx);
    break_points = find(diff_idx > 1);
    start_arcs = [thrust_idx(1), thrust_idx(break_points + 1)];
    end_arcs = [thrust_idx(break_points), thrust_idx(end)];
    
    for k = 1:length(start_arcs)
        idx_s = start_arcs(k);
        idx_e = end_arcs(k);
        idx_m = round((idx_s + idx_e)/2); % Punto medio
        
        % Sovrascrivi parte spinta con linea rossa
        plot(x_plot(idx_s:idx_e), y_plot(idx_s:idx_e), 'LineWidth', 2, 'Color',[0.4 0.4 0.4],'HandleVisibility', 'on');
        
        % Markers On/Off
        plot(x_plot(idx_s), y_plot(idx_s), '^','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor', [0 0.7 0], 'MarkerSize', 6, 'HandleVisibility', 'on');
        plot(x_plot(idx_e), y_plot(idx_e), 's','MarkerEdgeColor',[0.7 0 0], 'MarkerFaceColor',[0.7 0 0], 'MarkerSize', 6, 'HandleVisibility', 'on');
        
        % --- Frecce direzione spinta: ~10 per arco ---
        N_arrows = 10;
        idx_list = unique(round(linspace(idx_s, idx_e, min(N_arrows, idx_e-idx_s+1))));
        
        for ia = idx_list
            u2 = DELTA(ia);                 % u2: pitch, positivo clockwise, in [-pi,pi]
            if ~isnan(u2)
        
                % wrap coerente in [-pi,pi] anche senza toolbox
                u2 = atan2(sin(u2), cos(u2));
        
                % Angolo inerziale (CCW da +x):
                % tangenziale = theta + pi/2, e u2 positivo clockwise => sottraggo u2
                ang = x2(ia) + pi/2 - u2;
        
                dx = cos(ang) * arrow_len;
                dy = sin(ang) * arrow_len;
        
                % Per non riempire la legenda con 100 frecce:
                hv = 'off';
                if (k==1) && (ia==idx_list(1)), hv = 'on'; end
        
                quiver(x_plot(ia), y_plot(ia), dx, dy, 0, ...
                    'Color',[0.7 0 0.7], 'LineWidth', 1, 'MaxHeadSize', 0.5, ...
                    'HandleVisibility', hv);
            end
        end

    end
end
plot(x_plot(idx_e), y_plot(idx_e),'o','MarkerEdgeColor', [0 0.8 0.8], 'MarkerFaceColor',[0 0.9 0.9],'HandleVisibility', 'off');      
title('Optimized Trajectory', 'Color', 'k');
xlabel('x [km]', 'Color', 'k'); ylabel('y [km]', 'Color', 'k');
legend('Uranus', 'Umbriel Orbit', 'Start', 'Coasting','Thrust','On','Off','Direction Arrow', 'Location', 'bestoutside', 'TextColor', 'k');
%%
% --- PLOT 3: ANGOLO DI SPINTA ---
figure(5);
plot(t/60/60/24, rad2deg(DELTA), '-','LineWidth',2.5,'Color',[0 0 0]);hold on;
plot(t/60/60/24, rad2deg(DELTA), '-','LineWidth',1.5,'Color',[0.7 0.7 0.7]);
title('Pitch Angle (Clockwise)');
xlabel('Time [days]'); ylabel('Delta [deg]');
xlim([0,max(t/60/60/24)])
ylim([min(rad2deg(DELTA)) max(rad2deg(DELTA))])
grid on;

% --- PLOT 4: VELOCITÀ ---
figure(6);
subplot(3,1,1)
plot(t(1:idx_stop), x3, '-', 'LineWidth', 2,'Color',[0.8 0.3 0.5]);
title('Radial Velocity (vr)'); grid on;
xlim([0,T_fin_fin])
subplot(3,1,2)
plot(t(1:idx_stop), x4, '-', 'LineWidth', 2,'Color',[0.3 0.5 0]);
title('Tangential Velocity (vtheta)'); ylabel('Velocity [km/s]') ; grid on;
xlim([0,T_fin_fin])
subplot(3,1,3)
plot(t(1:idx_stop), sqrt(x3((1:idx_stop)).^2 + x4((1:idx_stop)).^2), 'k-', 'LineWidth', 2);
title('Velocity Module'); xlabel('Time [s]'); grid on;
xlim([0,T_fin_fin])


% --- PLOT 5: SWITCH FUNCTION E PROFILO DI SPINTA ---
figure(7);
hold on;
plot(t, CONTROL*1e3,'LineWidth', 2,'Color',[0 0.5 0]);
hold on;
plot(t, Switch_Function/10, '-', 'LineWidth', 2,'Color',[0.6 0.3 0.9]);
hold on; yline(0, 'k--', 'LineWidth', 2);
title('5. Switch Function and Control Profile');
xlabel('Time [s]');
ylabel('u_T [m/s^2]')
xlim([0,max(t)])
legend('Control Profile','Switch Function*1/10', 'Switch (0)');
grid on;

Active_time=CONTROL(:, any(CONTROL, 1));
Propulsion_Time=T_fin_fin*(length(Active_time)/length(CONTROL));
DeltaVFIN=(Propulsion_Time*gamma_max+norm(x3(end)));

fprintf('DeltaV: %2.4f km/s \n',DeltaVFIN)
fprintf('H(tf): %2.10e \n',H_star(end))

DV1 = sqrt(mu_ur/R_start) * (sqrt( 2*r_umbriel / (r_umbriel+R_start) ) - 1)
DV2 = sqrt(mu_ur/r_umbriel) * (1 - sqrt( 2*R_start / (r_umbriel+R_start) ) )
abs(DV1)+abs(DV2)