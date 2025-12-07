function plot_trajectory_TU(soluzione)
% PLOT_TRAJECTORY_TU
% Visualizza la traiettoria interplanetaria Terra -> Marte -> Giove -> Urano
% ricostruendo gli archi orbitali dai dati della struttura 'best_soluzione'.
%
% INPUT:
%   soluzione: la struct 'best_soluzione' estratta dalla cella di output
%              della funzione TERRA_MARTE_GIOVE_URANO.

    %% 1. Inizializzazione Grafica
    figure('Name', 'Missione Multi-Flyby: T-M-J-U', 'Color', 'w', 'NumberTitle', 'off');
    hold on;
    grid on;
    axis equal;
    xlabel('x_{ecliptic} [km]');
    ylabel('y_{ecliptic} [km]');
    title('Traiettoria Interplanetaria: Terra \rightarrow Urano');
    
    % Recupero GM del Sole (necessario per propagare le orbite per il plot)
    % Assumiamo che i kernel siano caricati nel MAIN.
    try
        gm_Sole = cspice_bodvrd('Sun', 'GM', 1);
    catch
        gm_Sole = 132712440041.939; % Valore di fallback se CSPICE fallisce
        warning('Impossibile recuperare GM_Sun da CSPICE. Uso valore di default.');
    end

    % Definizione Colori
    c_Sun = [1, 0.85, 0];
    c_Earth = [0, 0.4, 1];
    c_Mars = [1, 0.2, 0];
    c_Jup = [0.8, 0.5, 0.2];
    c_Ura = [0, 0.8, 0.8];
    c_Traj = [0.2, 0.2, 0.2]; % Grigio scuro per la traiettoria

    % Plot del Sole
    plot(0, 0, 'o', 'MarkerSize', 14, 'MarkerFaceColor', c_Sun, 'Color', c_Sun, 'DisplayName', 'Sole');

    %% 2. Plotting delle Tratte
    % La logica è: prendiamo lo stato vettore finale di ogni tratta (arrivo al pianeta)
    % e propaghiamo all'indietro per il tempo 'delta_t' per disegnare l'arco.

    % --- LEG 1: TERRA -> MARTE ---
    if isfield(soluzione, 'Earth_Mars_tr')
        tratta1 = soluzione.Earth_Mars_tr;
        % Stato all'arrivo (Marte)
        rf = tratta1.r'; % Assicuro vettore colonna
        vf = tratta1.v'; 
        dt = tratta1.delta_t;
        
        % Disegna arco e restituisce punto di partenza (Terra)
        r_start_T = plot_arc_backwards(rf, vf, dt, gm_Sole, c_Traj, 'Terra-Marte');
        
        % Plot Marker Pianeti
        plot_planet(r_start_T, c_Earth, 'Partenza Terra');
        plot_planet(rf, c_Mars, 'Flyby Marte');
    else
        warning('Campo Earth_Mars_tr mancante.');
    end

    % --- LEG 2: MARTE -> GIOVE ---
    if isfield(soluzione, 'Mars_Jupiter_tr')
        tratta2 = soluzione.Mars_Jupiter_tr;
        rf = tratta2.r';
        vf = tratta2.v';
        dt = tratta2.delta_t;
        
        plot_arc_backwards(rf, vf, dt, gm_Sole, c_Traj, 'Marte-Giove');
        plot_planet(rf, c_Jup, 'Flyby Giove');
    end

    % --- LEG 3: GIOVE -> URANO ---
    if isfield(soluzione, 'Jupiter_Uranus_tr')
        tratta3 = soluzione.Jupiter_Uranus_tr;
        rf = tratta3.r';
        vf = tratta3.v';
        dt = tratta3.delta_t;
        
        plot_arc_backwards(rf, vf, dt, gm_Sole, c_Traj, 'Giove-Urano');
        plot_planet(rf, c_Ura, 'Arrivo Urano');
    end

    %% 3. Riferimenti (Orbite dei pianeti approssimate)
    % Disegno cerchi tratteggiati per dare contesto visivo
    theta = linspace(0, 2*pi, 300);
    plot(1.496e8*cos(theta), 1.496e8*sin(theta), '--', 'Color', [c_Earth, 0.3], 'HandleVisibility', 'off');
    plot(2.279e8*cos(theta), 2.279e8*sin(theta), '--', 'Color', [c_Mars, 0.3], 'HandleVisibility', 'off');
    plot(7.785e8*cos(theta), 7.785e8*sin(theta), '--', 'Color', [c_Jup, 0.3], 'HandleVisibility', 'off');
    plot(2.871e9*cos(theta), 2.871e9*sin(theta), '--', 'Color', [c_Ura, 0.3], 'HandleVisibility', 'off');

    legend('show', 'Location', 'bestoutside');
    hold off;

end

%% --- FUNZIONI AUSILIARIE LOCALI ---

function r_start = plot_arc_backwards(r_final, v_final, dt, gm, col, name)
    % PLOT_ARC_BACKWARDS
    % Propaga un'orbita a due corpi all'indietro dal punto finale
    % Ritorna la posizione al punto iniziale.
    
    n_points = 100;
    times = linspace(0, -dt, n_points); % Tempi negativi (indietro)
    
    vals = zeros(6, n_points);
    
    % Se cspice è disponibile usiamo prop2b per massima precisione
    % Altrimenti fallback su propagatore kepleriano semplice
    try
        % cspice_prop2b propaga lo stato vettore (state) per un tempo (t)
        % state è [6x1], t può essere array
        state0 = [r_final; v_final];
        for k = 1:n_points
            vals(:,k) = cspice_prop2b(gm, state0, times(k));
        end
    catch
        error('CSPICE non trovato o cspice_prop2b fallito. Assicurati di caricare Mice.');
    end
    
    % Estrazione coordinate
    X = vals(1,:);
    Y = vals(2,:);
    Z = vals(3,:);
    
    % Plot Traiettoria
    plot3(X, Y, Z, 'Color', col, 'LineWidth', 1.5, 'DisplayName', name);
    
    % Ritorna il punto iniziale (l'ultimo calcolato poichè andiamo indietro)
    r_start = vals(1:3, end);
end

function plot_planet(r, col, name)
    % Helper per disegnare il pallino del pianeta
    plot3(r(1), r(2), r(3), 'o', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', col, ...
        'Color', 'k', ...
        'DisplayName', name);
end