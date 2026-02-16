%% OTTIMIZZAZIONE TRAIETTORIA SATURNO - ENCELADO
clc; clear variables; close all;

%% 1. DATI GLOBALI E FISICI
global gamma_max W R_start r_umbriel mu_ur g0



% --- Parametri Fisici ---
mu_ur = 5793939;                 % [km^3/s^2] Parametro gravitazionale Saturno
r_ur = 25362;                     % [km] Raggio Saturno
r_umbriel = 266000;             % [km] Raggio orbita Encelado
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
% --- Parametri Motore ---
% Spinta specifica (accelerazione)
% gamma_max = 0.8 / (9.81 * 1300); % [km/s^2] (approx 6.27e-5 m/s^2)
% W = gamma_max * 1300;            % Peso fittizio per equazione massa


% --- Definizione Orbita di Partenza (ESTERNA) ---
% Partiamo da un'orbita circolare più alta di quella di Umbriel
R_start = r_umbriel + 400000;   % Partenza a rp~416000 km

%% 2. CONDIZIONI AL CONTORNO (SCALING)
% Stato Iniziale (Start)
r0 = R_start;

e0=0.6;
a0=r0/(1-e0);

theta0 = 0; 
vr0 = 0;
vtheta0 = sqrt((mu_ur*(1+e0))/r0); 

% Stato Finale Target (Encelado) - Usato solo per riferimento qui
rf_target = r_umbriel;

% normalizzazione
L  = R_start;
V  = sqrt(mu_ur*(1+0.6) / L);
T  = L / V;
A  = V / T;  

%% 3. IMPOSTAZIONE OTTIMIZZAZIONE (FIRST GUESS & BOUNDS)


% Stima del Tempo di Volo ( vers precedente !!!)
% Stima conservativa: circa 70 giorni per scendere dolcemente
t_fin_guess = 3669600;    
% Cf_guess = gamma_max * t_fin_guess * 0.5; 

% Stima del Tempo di Volo
% Stima basata su circa 1.5 anni
% t_fin_guess = 5.0e7;    


% --- BOUNDS (LIMITI) AGGIORNATI ---
% Allargo i limiti su Pvr (secondo e terzo elemento) per permettere
% manovre più aggressive per evitare il pianeta.
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
  'MaxIterations', 30, ...
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

% --- LOGICA DI STOP (TAGLIO DATI SENZA EVENTO) ---
% Cerco il primo indice in cui il raggio scende sotto o è uguale al raggio di Encelado.
% Siccome stiamo scendendo dall'alto, questo è il punto di contatto.
%idx_stop = find(integrator_output(:,1) <= r_umbriel, 1, 'first'); 
idx_stop = size(integrator_output(:,1),1);
t = t * T;
day_stop=floor(t(idx_stop)/60/60/24);
hour_stop=floor((t(idx_stop)/60/60)-day_stop*24);
T_fin_fin=t(idx_stop);

if ~isempty(idx_stop)
    fprintf('Target Raggiunto a t =%d day ,%d hour (su un totale stimato di %.2f s).\n',day_stop,hour_stop,t_fin);
    fprintf('Troncamento della traiettoria in eccesso.\n');
    
    % Taglio tutti i vettori a questo indice
    t = t(1:idx_stop);
    integrator_output = integrator_output(1:idx_stop, :);
else
    warning('La traiettoria non ha toccato l''orbita di Umbriel nel tempo previsto.');
end
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
    Switch_Function = [Switch_Function; SF*100];
    
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
title('1. Hamiltonian (Constant)');
xlabel('Time [s]'); ylabel('H');
xlim([0,max(t)])
grid on;
%%
% --- PLOT 2: TRAIETTORIA (Advanced) ---
figure(2); hold on; grid on; axis equal;
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
arrow_len = 40000; % Lunghezza frecce (regolabile)

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
        
        % Frecce Direzione Spinta (Start, Mid, End)
        indices_arrow = [idx_s, idx_m, idx_e];
        for ia = indices_arrow
             if ~isnan(DELTA(ia))
                ang_inertial = x2(ia) + DELTA(ia); % Angolo assoluto
                u = sin(ang_inertial); v = cos(ang_inertial);
                quiver(x_plot(ia), y_plot(ia), u*arrow_len, v*arrow_len, 0, ...
                       'Color',[0.7 0 0.7], 'LineWidth', 1, 'MaxHeadSize', 0.5, 'HandleVisibility', 'on');
             end
        end
    end
end
plot(x_plot(idx_e), y_plot(idx_e),'o','MarkerEdgeColor', [0 0.8 0.8], 'MarkerFaceColor',[0 0.9 0.9],'HandleVisibility', 'off');      
title('2. Optimized Trajectory', 'Color', 'k');
xlabel('x [km]', 'Color', 'k'); ylabel('y [km]', 'Color', 'k');
legend('Uranus', 'Umbriel Orbit', 'Start', 'Coasting','Thrust','On','Off','Direction Arrow', 'Location', 'bestoutside', 'TextColor', 'k');
%%
% --- PLOT 3: ANGOLO DI SPINTA ---
figure(3);
plot(t/60/60/24, rad2deg(DELTA), '-','LineWidth',2.5,'Color',[0 0 0]);hold on;
plot(t/60/60/24, rad2deg(DELTA), '-','LineWidth',1.5,'Color',[0.7 0.7 0.7]);
title('3. Thrust Angle (w.r.t the local)');
xlabel('Time [days]'); ylabel('Delta [deg]');
xlim([0,max(t/60/60/24)])
ylim([min(rad2deg(DELTA)) max(rad2deg(DELTA))])
grid on;

% --- PLOT 4: VELOCITÀ ---
figure(4);
subplot(3,1,1)
plot(t(1:idx_stop), x3, '-', 'LineWidth', 2,'Color',[0.8 0.3 0.5]);
title('Radial Velocity (vr)'); grid on;
xlim([0,T_fin_fin])
subplot(3,1,2)
plot(t(1:idx_stop), x4, '-', 'LineWidth', 2,'Color',[0.3 0.5 0]);
title('Tangential Velocity (vtheta)'); grid on;
xlim([0,T_fin_fin])
subplot(3,1,3)
plot(t(1:idx_stop), sqrt(x3((1:idx_stop)).^2 + x4((1:idx_stop)).^2), 'k-', 'LineWidth', 2);
title('Velocity Module'); xlabel('Time [s]'); grid on;
xlim([0,T_fin_fin])


% --- PLOT 5: SWITCH FUNCTION E PROFILO DI SPINTA ---
figure(5);
hold on;
plot(t, CONTROL*1e6,'LineWidth', 2,'Color',[0 0.5 0]);
hold on;
plot(t, Switch_Function, '-', 'LineWidth', 2,'Color',[0.6 0.3 0.9]);
hold on; yline(0, 'k--', 'LineWidth', 2);
title('5. Switch Function and Control Profile');
xlabel('Time [s]');
xlim([0,max(t)])
legend('Control Profile','100*Switch Function', 'Switch (0)');
grid on;

Active_time=CONTROL(:, any(CONTROL, 1));
Propulsion_Time=T_fin_fin*(length(Active_time)/length(CONTROL));
DeltaVFIN=(Propulsion_Time*gamma_max+norm(x3(end)))*1000;

fprintf('DeltaV: %2.4f m/s \n',DeltaVFIN)
fprintf('H(tf): %2.10e \n',H_star(end))
