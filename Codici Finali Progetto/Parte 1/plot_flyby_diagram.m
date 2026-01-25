function gravity_assist_diagram()
    % Setup della figura
    figure('Color', 'w', 'Position', [100, 100, 800, 500]);
    hold on;
    axis equal;
    box on;
    grid off;
    
    % --- DEFINIZIONE PARAMETRI GEOMETRICI ---
    % Coordinate e magnitudini stimate dall'immagine
    origin = [0, 0];          % Il "pianeta" dove i vettori v_inf divergono
    Vp_start = [8, 0];        % La coda del vettore Vp (arancione)
    v_inf_mag = 4.5;          % Magnitudine v_infinito
    
    % Angoli (rispetto all'asse orizzontale sinistro, alpha=0)
    % Nota: L'immagine mostra alpha=0 a sinistra (180 gradi matematici standard)
    % Quindi dobbiamo lavorare con questo sistema di riferimento.
    angle_minus_deg = 55;     % Angolo del vettore marrone
    angle_plus_deg = 135;     % Angolo del vettore blu
    
    % Calcolo coordinate punte dei vettori
    % Convertiamo gli angoli in radianti e coordinate cartesiane
    % (Invertiamo la x perché l'origine degli angoli è a sinistra)
    vec_inf_minus = [-v_inf_mag * cosd(angle_minus_deg), v_inf_mag * sind(angle_minus_deg)];
    vec_inf_plus  = [-v_inf_mag * cosd(angle_plus_deg),  v_inf_mag * sind(angle_plus_deg)];
    
    % --- DEFINIZIONE COLORI (RGB) ---
    col_orange = [0.93, 0.69, 0.13]; % Vp
    col_blue   = [0.0, 0.3, 0.6];    % v_inf+
    col_brown  = [0.8, 0.6, 0.4];    % v_inf-
    col_red    = [0.8, 0.1, 0.1];    % V+, V- (Heliocentric)
    col_green  = [0.4, 0.8, 0.2];    % Angoli
    col_dark_green = [0.0, 0.5, 0.3]; % Delta
    
    % --- DISEGNO VETTORI ---

     % --- DISEGNO ASSI E RIFERIMENTI ---
    % Linea grigia di base
    line([-5, 10], [0, 0], 'Color', [0.8 0.8 0.8], 'LineWidth', 2);
    
    % 1. Vettore Vp (Arancione) - Da destra verso l'origine
    % Disegniamo una linea spessa arancione da Vp_start a origin
    draw_vector2(Vp_start, origin - Vp_start, col_orange, '$\bar{V}_p$');
    
    % 2. Vettori v_infinito (Dall'origine verso l'esterno)
    % Marrone (v_inf -)
    draw_vector(origin, vec_inf_minus, col_brown, '$\bar{v}_{\infty}^+$');
    % Blu (v_inf +)
    draw_vector(origin, vec_inf_plus, col_blue, '$\bar{v}_{\infty}^-$');
    
    % 3. Vettori Eliocentrici (Rossi)
    % Collegano la coda di Vp alle punte di v_inf
    % V- (Rosso verso marrone)
    target_minus = origin + vec_inf_minus;
    draw_vector2(Vp_start, target_minus - Vp_start, col_red, '$\bar{V}^+$');
    
    % V+ (Rosso verso blu)
    target_plus = origin + vec_inf_plus;
    draw_vector2(Vp_start, target_plus - Vp_start, col_red, '$\bar{V}^-$');
    
    
    % Punti neri su alpha=0 e alpha=pi
    plot(-4, 0, 'k.', 'MarkerSize', 20);
    plot(Vp_start(1), 0, 'k.', 'MarkerSize', 15); % Punto su alpha=pi approssimativo
    
    text(-4.2, -0.5, '$\alpha=0$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
    text(Vp_start(1), -0.5, '$\alpha=\pi$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
    
    % --- DISEGNO ARCHI ANGOLI ---
    
    % Arco alpha- (Verde chiaro) - Da asse sx a vettore marrone
    draw_arc(origin, 2.5, 180, 180-angle_minus_deg, col_green);
    % Etichetta alpha-
    text(-3, 1.5, '$\alpha^+$', 'Interpreter', 'latex', 'Color', col_green, 'FontSize', 16);
    
    % Arco alpha+ (Verde chiaro) - Da asse sx a vettore blu
    draw_arc(origin, 3.5, 180, 180-angle_plus_deg+2, col_green);
    % Etichetta alpha+
    text(-3.8, 2.2, '$\alpha^-$', 'Interpreter', 'latex', 'Color', col_green, 'FontSize', 16);
    
    % Arco Delta (Verde scuro) - Tra vettore blu e marrone
    draw_arc(origin, 1.5, 180-angle_minus_deg, 180-angle_plus_deg, col_dark_green);
    % Etichetta Delta
    text(-0.5, 1.0, '$\delta$', 'Interpreter', 'latex', 'Color', col_dark_green, 'FontSize', 18);

    % --- PULIZIA GRAFICO ---
    xlim([-6, 10]);
    ylim([-1, 6]);
    set(gca, 'XTickLabel', [], 'YTickLabel', []); % Nascondi numeri assi
end

function draw_vector(start_pt, components, color_vec, label_str)
    % Wrapper per disegnare vettori con quiver e gestire le etichette
    q = quiver(start_pt(1), start_pt(2), components(1), components(2), ...
        0, 'Color', color_vec, 'LineWidth', 3, 'MaxHeadSize', 0.15);
    
    % Calcolo posizione label (a metà vettore o alla punta con un offset)
    % Qui posizioniamo vicino alla punta
    tip = start_pt + components;
    
    % Logica offset manuale per assomigliare all'immagine
    offset = [0, 0.3]; 
    if contains(label_str, 'V_p')
        pos = start_pt + components/2; % Vp etichetta al centro
        offset = [0, 0.3];
    elseif contains(label_str, 'V^+')
        pos = tip; 
        offset = [-0.5, 0.3]; % Sposta un po' a sinistra
    elseif contains(label_str, 'v_{\infty}^+')
        pos = tip; 
        offset = [0.1, 0.3];
    else
        pos = tip;
    end
    
    text(pos(1)+offset(1), pos(2)+offset(2), label_str, ...
        'Interpreter', 'latex', 'Color', color_vec, 'FontSize', 16, 'FontWeight', 'bold');
end

function draw_vector2(start_pt, components, color_vec, label_str)
    % Wrapper per disegnare vettori con quiver e gestire le etichette
    q = quiver(start_pt(1), start_pt(2), components(1), components(2), ...
        0, 'Color', color_vec, 'LineWidth', 3, 'MaxHeadSize', 0.15);
    
    % Calcolo posizione label (a metà vettore o alla punta con un offset)
    % Qui posizioniamo vicino alla punta
    tip = start_pt + components;
    
    % Logica offset manuale per assomigliare all'immagine
    offset = [-0.3, -0.6]; 
    if contains(label_str, 'V_p')
        pos = start_pt + components/2; % Vp etichetta al centro
        offset = [1, 1];
    elseif contains(label_str, 'V^+')
        pos = tip; 
        offset = [-0.5, -0.3]; % Sposta un po' a sinistra
    elseif contains(label_str, 'v_{\infty}^+')
        pos = tip; 
        offset = [0.1, 0.3];
    else
        pos = tip;
    end
    
    text(pos(1)+offset(1), pos(2)+offset(2), label_str, ...
        'Interpreter', 'latex', 'Color', color_vec, 'FontSize', 16, 'FontWeight', 'bold');
end

function draw_arc(center, radius, start_angle_deg, end_angle_deg, color_vec)
    % Disegna un arco e aggiunge le frecce alle estremità
    t = linspace(start_angle_deg, end_angle_deg, 50);
    x = center(1) + radius * cosd(t);
    y = center(2) + radius * sind(t);
    plot(x, y, 'Color', color_vec, 'LineWidth', 2);
    
    % Aggiungi freccia alla fine dell'arco
    % Calcoliamo la direzione tangente finale per orientare la freccia
    % Metodo semplice: usiamo quiver sull'ultimo punto
    u_arrow = x(end) - x(end-1);
    v_arrow = y(end) - y(end-1);
    quiver(x(end), y(end), u_arrow, v_arrow, 0.5, ...
        'Color', color_vec, 'LineWidth', 2, 'MaxHeadSize', 0.5, 'ShowArrowHead', 'on');
end