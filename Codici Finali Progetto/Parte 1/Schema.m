% Pulizia dell'ambiente
clear; clc; close all;

% --- CONFIGURAZIONE DATI ---

% Distanze medie dal Sole in UA (Unità Astronomiche)
r_terra = 1.0;
r_marte = 1.52;
r_giove = 5.2;
r_urano = 19.2;

% Vettori dei raggi e degli angoli (Theta) per i fly-by
% Gli angoli simulano la posizione progressiva durante il viaggio
r_punti = [r_terra, r_marte, r_giove, r_urano];
theta_punti = [0, 1.2, 3.5, 5.8]; % Radianti

% Interpolazione (Spline) per creare una curva morbida della traiettoria
theta_fine = linspace(min(theta_punti), max(theta_punti), 500);
r_fine = spline(theta_punti, r_punti, theta_fine);

% Conversione in coordinate Cartesiane per la traiettoria
x_traj = r_fine .* cos(theta_fine);
y_traj = r_fine .* sin(theta_fine);

% Coordinate esatte dei pianeti (punti di incontro)
x_planets = r_punti .* cos(theta_punti);
y_planets = r_punti .* sin(theta_punti);

% --- CREAZIONE DEL GRAFICO ---

% Creazione figura con sfondo nero
figure('Color', 'w', 'Name', 'Traiettoria Interplanetaria');
hold on;
axis equal; % Mantiene le proporzioni corrette (non deforma i cerchi)

% 1. Disegno le orbite di riferimento (cerchi grigi tenui)
angoli_orbita = linspace(0, 2*pi, 300);
raggi_pianeti = [r_terra, r_marte, r_giove, r_urano];

for r = raggi_pianeti
    x_orb = r * cos(angoli_orbita);
    y_orb = r * sin(angoli_orbita);
    plot(x_orb, y_orb, '--', 'Color', [0 0 0], 'LineWidth', 0.5);
end

% 2. Disegno la Traiettoria della Sonda (Bianca)
plot(x_traj, y_traj, 'k-', 'LineWidth', 2);

% 3. Disegno il Sole (Giallo Scuro)
% Colore RGB: [0.7, 0.5, 0] -> Ocra/Giallo scuro
scatter(0, 0, 400, 'MarkerFaceColor', [0.7, 0.5, 0], ...
    'MarkerEdgeColor', [1, 0.6, 0], 'LineWidth', 1.5);

% 4. Disegno i Pianeti (con colori specifici RGB)

% Terra (Azzurra) - Start
scatter(x_planets(1), y_planets(1), 80, 'filled', ...
    'MarkerFaceColor', [0, 0.7, 1]); % Deep Sky Blue

% Marte (Rosso) - Flyby 1
scatter(x_planets(2), y_planets(2), 60, 'filled', ...
    'MarkerFaceColor', [1, 0.2, 0.2]); % Rosso

% Giove (Marrone Chiaro) - Flyby 2
scatter(x_planets(3), y_planets(3), 200, 'filled', ...
    'MarkerFaceColor', [0.8, 0.6, 0.4]); % Peru/Marrone chiaro

% Urano (Celeste Spento) - Arrivo
scatter(x_planets(4), y_planets(4), 120, 'filled', ...
    'MarkerFaceColor', [0.7, 0.9, 0.95]); % Powder Blue/Celeste pallido

% --- ESTETICA E TESTI ---

% Aggiunta etichette (Text)
text(x_planets(1)+0.5, y_planets(1)+0.5, 'Start (Earth)', 'Color', 'k', 'FontSize', 9);
text(x_planets(2)+0.5, y_planets(2), 'Gravity Assist 1', 'Color', [1, 0.5, 0.5], 'FontSize', 8);
text(x_planets(3)+1, y_planets(3), 'Gravity Assist 2', 'Color', [0.9, 0.8, 0.6], 'FontSize', 8);
text(x_planets(4)+1, y_planets(4), 'Target (Uranus)', 'Color', [0, 0, 0], 'FontSize', 9);

% Impostazioni Assi (Sfondo nero e testo grigio)
ax = gca;
ax.Color = 'w';           % Sfondo assi nero
%grid on;
ax.XColor = [0 0 0]; % Colore numeri assi X
ax.YColor = [0 0 0]; % Colore numeri assi Y
xlabel('X - (UA)', 'Color', [0 0 0]);
ylabel('Y - (UA)', 'Color', [0 0 0]);
title('Trajectory: Earth -> Mars -> Jupiter -> Uranus', 'Color', 'w', 'FontSize', 12);

% Zoom per inquadrare tutto
xlim([-8 22]);
ylim([-15 5]);

hold off;