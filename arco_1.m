function [V_inf_terra,V_ext_terra_SoI,True_anomaly_Mars_arrival,V_minus,Delta_i_ing] = arco_1(delta_V,h_LEO)
%INPUT:1°DeltaV 2°h_LEO
%OUTPUT:1°V inifinito SC usciti dalla SoI della terra.2°V Assoluta SC
%usciti dalla SoI terra. 3°Anomalia vera all'arrivo su marte. 4°V SC
%Assoluta all'incontro con SoI di marte. 5°DELTA iniziale.


%% Definition of the problem constants (all in km and km/s)
MuM = 4.282837e4;              % Mars gravitational parameter [km^3/s^2]
MuE = 3.986004418e5;           % Earth gravitational parameter [km^3/s^2]
MuS = 1.32712440018e11;        % Sun gravitational parameter [km^3/s^2]

Re  = 1.495978707e8;           % Earth's orbital semimajor axis [km]
Rm  = 2.2794382e8;             % Mars' orbital semimajor axis [km]
R_EM = Rm - Re;                % Average Earth–Mars distance [km]

ve = sqrt(MuS/Re);             % Earth's orbital speed around the Sun [km/s]
vm = sqrt(MuS/Rm);             % Mars' orbital speed around the Sun [km/s]

%% Orbite coplanari (per rappresentazione)
theta = linspace(0, 2*pi, 500);   % angolo (rad)

xE = Re * cos(theta);
yE = Re * sin(theta);

xM = Rm * cos(theta);
yM = Rm * sin(theta);


%% Input: condizione di partenza da LEO

r_leo = Re + h_LEO;                % [km]
v_c = sqrt(MuE / r_leo);           % velocità circolare in LEO
v_p = v_c + delta_V;               % velocità dopo il delta-V
v_esc = sqrt(2 * MuE / r_leo);     % velocità di fuga

% Calcolo v_inf
if v_p < v_esc
    v_inf = 0;
else
    v_inf = sqrt(v_p^2 - v_esc^2);
end
V_inf_terra=v_inf;
V_ext_terra_SoI=v_inf+ve;

%% Parametri dell'arco di trasferimento
v_minus = v_inf + ve;                         % velocità eliocentrica al di fuori della SOI terrestre
a = (2/Re - v_minus^2 / MuS)^(-1);            % semiasse maggiore del trasferimento
e = 1 - Re / a;                               % eccentricità

% Controlli
if e < 0
    error('ERROR: eccentricità negativa. Controlla la v_inf o il delta-V.');
end

value = (a*(1 - e^2)/(e*Rm)) - 1/e;
if value < -1 || value > 1
    error('ERROR: valore di true anomaly complesso. Controlla la v_inf o il delta-V.');
end

theta_star_plus = acos(value);  % true anomaly all'arrivo
True_anomaly_Mars_arrival= rad2deg(theta_star_plus);

%% Coordinate della traiettoria di trasferimento
theta_arc = linspace(0, theta_star_plus, 1000);
r_arc = (a*(1 - e^2)) ./ (1 + e*cos(theta_arc));
x_arc = r_arc .* cos(theta_arc);
y_arc = r_arc .* sin(theta_arc);

%% Calcolo della velocità di arrivo su Marte
% momento angolare specifico dell'orbita eliocentrica di trasferimento
h = sqrt(MuS * a * (1 - e^2));

% componenti della velocità del S/C (eliocentrica) al true anomaly theta_star_plus
v_r = (MuS / h) * e * sin(theta_star_plus);           % componente radiale [km/s]
v_theta = (MuS / h) * (1 + e * cos(theta_star_plus));% componente tangenziale [km/s]

v_sc_arr = sqrt(v_r^2 + v_theta^2);                   % velocità eliocentrica dello S/C all'arrivo [km/s]

% velocità eliocentrica di Marte (tangenziale, in modulo)
% (già calcolata come vm sopra) - direzione tangenziale positiva
% velocità relativa (iperbolic excess) rispetto a Marte
v_rel_r = v_r;                 % Marte ha componente radiale nulla nell'approssimazione circolare coplanare
v_rel_theta = v_theta - vm;
v_inf_Mars = sqrt(v_rel_r^2 + v_rel_theta^2);

% Stampa risultati (condizioni di arrivo su marte)
V_minus= v_sc_arr;

% angolo con segno: positivo => v_sc ha componente radiale uscente
delta_i = atan2( v_r, (-v_rel_theta) );   % rad
delta_i_deg = rad2deg(delta_i);
Delta_i_ing=delta_i_deg;

%% Grafico finale combinato
% figure('Color', 'w'); hold on; grid on; axis equal;
% 
% % Sole
% plot(0, 0, 'yo', 'MarkerFaceColor','y', 'MarkerSize', 10);
% 
% % Orbite Terra e Marte
% plot(xE, yE, 'b-', 'LineWidth', 1.2);
% plot(xM, yM, 'r--', 'LineWidth', 1.2);
% 
% % Arco di trasferimento
% plot(x_arc, y_arc, 'k', 'LineWidth', 2);
% 
% % Pianeti
% plot(Re, 0, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 6);
% plot(Rm*cos(theta_star_plus), Rm*sin(theta_star_plus), 'ro', ...
%     'MarkerFaceColor','r', 'MarkerSize', 6);
% 
% title('Trasferimento interplanetario Terra → Marte (frame eliocentrico)');
% xlabel('x [km]');
% ylabel('y [km]');
% legend('Sole','Orbita Terra','Orbita Marte','Arco di trasferimento','Terra','Marte', ...
%        'Location','bestoutside');
% 
% axis([-2.5e8 2.5e8 -2.5e8 2.5e8]);
end

