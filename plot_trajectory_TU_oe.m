function plot_trajectory_TU_oe(soluzione)
% PLOT_TRAJECTORY_TU_OE
% Visualizza la traiettoria usando l'equazione polare della conica.
% INCLUDE: Plot della Terra alla partenza.

    %% 1. Setup Grafico
    figure('Name', 'Traiettoria (Parametri Orbitali)', 'Color', 'w', 'NumberTitle', 'off');
    hold on; grid on; axis equal;
    xlabel('X Ecliptic [km]'); ylabel('Y Ecliptic [km]');
    title('Traiettoria Multi-Flyby: Terra \rightarrow Urano');
    
    % --- Definizione Colori ---
    c_Sun   = [1, 0.85, 0];    % Giallo Sole
    c_Earth = [0, 0, 1];       % Blu Terra
    c_Mars  = [1, 0, 0];       % Rosso Marte
    c_Jup   = [0.8, 0.5, 0.2]; % Arancione/Marrone Giove
    c_Ura   = [0, 0.8, 0.8];   % Ciano Urano
    c_Traj  = [0.2, 0.2, 0.2]; % Grigio scuro Traiettoria
    
    % Disegna il Sole al centro
    plot(0, 0, 'o', 'MarkerSize', 14, 'MarkerFaceColor', c_Sun, 'Color', c_Sun, 'DisplayName', 'Sole');

    %% 2. Plot delle Tratte e dei Pianeti
    
    % --- TRATTA 1: TERRA -> MARTE ---
        % Disegna l'arco e RECUPERA il punto iniziale (Terra)
        r_start_T = draw_segment(soluzione.Earth_Mars_tr, c_Traj, 'Terra-Marte');
       
    % --- TRATTA 2: MARTE -> GIOVE ---
        draw_segment(soluzione.Mars_Jupiter_tr, c_Traj, 'Marte-Giove');
        
    % --- TRATTA 3: GIOVE -> URANO ---
        draw_segment(soluzione.Jupiter_Uranus_tr, c_Traj, 'Giove-Urano');  
        
    % --- PLOT PIANETI ---
        % Terra
        plot_planet(r_start_T, c_Earth, 'Terra');
        
        % Marte (Arrivo tratta 1)
        plot_planet(soluzione.Earth_Mars_tr.r, c_Mars, 'Marte');

        % Giove (Arrivo tratta 2)
        plot_planet(soluzione.Mars_Jupiter_tr.r, c_Jup, 'Giove');
    
        % Urano (Arrivo tratta 3)
        plot_planet(soluzione.Jupiter_Uranus_tr.r, c_Ura, 'Urano');
    
    %% 3. Orbite di Riferimento
    
    plot_ref_orbit(1.496e8, c_Earth); % Terra orbit
    plot_ref_orbit(2.279e8, c_Mars);  % Marte orbit
    plot_ref_orbit(7.785e8, c_Jup);   % Giove orbit
    plot_ref_orbit(2.871e9, c_Ura);   % Urano orbit

    legend('show', 'Location', 'bestoutside');
    hold off;
end

%% --- FUNZIONI AUSILIARIE ---

function r_start_vec = draw_segment(struct_tratta, color, name)
    % Disegna l'arco di conica e restituisce il vettore posizione iniziale [x;y;0]
    
    info = struct_tratta.info;
    r_final_vec = struct_tratta.r; % Serve per calcolare l'orientamento
    
    % 1. Recupero Parametri
    a = info.a;
    e = info.e;
    th1 = info.theta1; % Anomalia vera partenza [rad]
    th2 = info.theta2; % Anomalia vera arrivo [rad]
    
    % Semilato retto
    p = a * (1 - e^2);
    
    % 2. Calcolo Orientamento (Argomento del pericentro 'w')
    % Logica: L'angolo polare finale geometrico è (w + th2).
    % Quindi w = (Angolo vettore r_finale) - th2
    angle_final_pos = atan2(r_final_vec(2), r_final_vec(1));
    w = angle_final_pos - th2; 
    
    % 3. Generazione Punti per il plot
    % Gestione wrap-around se th2 < th1 (es. passa per perielio)
    if th2 < th1
         nu = linspace(th1, th2, 100); % Assumiamo continuità matematica o la gestiamo linearmente
    else
         nu = linspace(th1, th2, 100);
    end
    
    % Equazione polare
    r = p ./ (1 + e .* cos(nu));
    
    % Coordinate nel piano perifocale
    x_perif = r .* cos(nu);
    y_perif = r .* sin(nu);
    
    % 4. Rotazione nel piano inerziale (Eclittica)
    % Applicazione rotazione di 'w' attorno a Z
    X_plot = x_perif * cos(w) - y_perif * sin(w);
    Y_plot = x_perif * sin(w) + y_perif * cos(w);
    
    % 5. Plot della linea
    plot(X_plot, Y_plot, 'LineWidth', 1.5, 'Color', color, 'DisplayName', name);
    
    % 6. Restituzione del punto iniziale (il primo punto calcolato)
    % Utile per plottare il pianeta di partenza (es. Terra)
    r_start_vec = [X_plot(1); Y_plot(1); 0];
end

function plot_planet(r, col, name)
    % Disegna un marker colorato nella posizione r
    plot(r(1), r(2), 'o', 'MarkerSize', 7, 'MarkerFaceColor', col, 'Color', 'k', 'DisplayName', name);
end

function plot_ref_orbit(radius, col)
    % Disegna orbita circolare di riferimento (tratteggiata)
    th = linspace(0, 2*pi, 300);
    % Aggiunge trasparenza al colore (quarto elemento alpha)
    col_transparent = [col, 0.3]; 
    if length(col) > 3, col_transparent = col; end % Se ha già alpha, lascialo
    
    plot(radius*cos(th), radius*sin(th), '--', 'Color', col_transparent, 'HandleVisibility', 'off');
end