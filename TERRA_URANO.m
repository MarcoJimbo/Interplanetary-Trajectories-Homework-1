function[SOLUZIONE] = TERRA_URANO(d_Terra,d_Urano,gm_Terra,gm_Urano,gm_Sole,SOI_Terra,SOI_Urano,v_Terra,v_Urano,rp_Terra,rp_Urano,E_Terra_0,E_Urano_0,ET_t0,k1,k12,dV_step)
%%% funzione che studia performance trasferta diretta TERRA URANO al variare della
%%% variabile di progetto dV_Terra (impulso erogato da orbita LEO).

%%% INPUT:
%%% d_Terra:         SCALAR [1x1] = distanza Terra-Sole [km] 
%%% d_Urano:         SCALAR [1x1] = distanza Urano-Sole [km]
%%% gm_Terra:        SCALAR [1x1] = parametro gravitazionale Terra [km^3/s^2]
%%% gm_Urano:        SCALAR [1x1] = parametro gravitazionale Urano [km^3/s^2]
%%% gm_Sole:         SCALAR [1x1] = parametro gravitazionale Sole [km^3/s^2]
%%% SOI_Terra:       SCALAR [1x1] = raggio SOI Terra [km]
%%% SOI_Urano:       SCALAR [1x1] = raggio SOI Urano [km]
%%% v_Terra:         SCALAR [1x1] = velocità Terra [km/s]
%%% v_Urano:         SCALAR [1x1] = velocità Urano [km/s]
%%% rp_Terra:        SCALAR [1x1] = raggio orbita LEO [km]
%%% rp_Urano:        SCALAR [1x1] = raggio orbita circolare target [km]
%%% E_Terra_0:       SCALAR [1x1] = anomalia eccentrica Terra a t0 [rad]
%%% E_Urano_0:       SCALAR [1x1] = anomalia eccentrica Urano a t0 [rad]
%%% ET_t0:           SCALAR [1x1] = t0 espresso in ET [s]
%%% k1:              SCALAR [1x1] = costante moltiplicativa che definisce dV_max = k1 * dV_min
%%% k12:             SCALAR [1x1] = dV_max se con dV=0 si raggiunge prossimo pianeta di missione [km/s]
%%% dV_step:         SCALAR [1x1] = risoluzione di variabile di progetto dV [km/s]

%%% OUTPUT:
%%% SOLUZIONE
%%%          .
%%%          .
%%%          .
%%%          .Earth_Exit
%%%                     .
%%%                     .
%%%                     .
%%%                     .
%%%                     .
%%%                     .info
%%%                          .
%%%                          .
%%%                          .
%%%                          .
%%%                          .
%%%                          .
%%%                          .
%%%                          .
%%%                          .
%%%                          .
%%%          .Int_Transfer
%%%                       .
%%%                       .
%%%                       .
%%%                       .
%%%                       .
%%%                       .
%%%                       .info
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%                            .
%%%          .Uranus_Entrance
%%%                          .
%%%                          .
%%%                          .
%%%                          .
%%%                          .info
%%%                               .
%%%                               .
%%%                               .
%%%                               .
%%%                               .
%%%                               .
%%%                               .


Vp_Terra = sqrt( gm_Terra / rp_Terra ); % [km/s]
[dV_Terra] = dV_setter(d_Terra,d_Urano,gm_Sole,gm_Terra,rp_Terra,Vp_Terra,k1,k12,dV_step);

% inizializzazione variabili di costo 
dV_tot = NaN(1,length(dV_Terra));
t_durata = zeros(size(dV_tot));
UTC_departure = strings(size(dV_tot));
UTC_arrival = strings(size(dV_tot));

% determinazione velocità di eccesso iperbolica in entrata SOI Urano [km/s]
% inizializzazione
SOLUZIONE = cell(1,length(dV_tot)); % cella contenente tutte le caratteristiche della missione 
v_inf = zeros(size(dV_tot)); % [km/s] velocità di eccesso iperbolico in entrata SOI Urano
r_SC_f = cell(1,length(dV_tot));
v_SC_f = cell(1,length(dV_tot));
X_SC_f = cell(1,length(dV_tot)); % statovettore dello Spacecraft all'encounter + informazioni su traiettoria interplanetaria
Earth_exit = cell(1,length(dV_tot)); % informazioni su traiettoria fuga da Terra 
Uranus_approach = cell(1,length(dV_tot)); % informazioni su traiettoria iperbolica in SOI Urano 
V_E = [ 0 , v_Terra , 0];
% calcolo
% poichè scelta sdr è arbitraria in questa fase scelgo sdr {O} centrato in Sole
% con x in direzione Terra e z in direzione perpendicolare al piano
% orbitale dei pianeti (ipotesi complanarità)
r_SC_i = d_Terra * [ 1 , 0 , 0 ]; % [km] posizione SC dopo uscita SOI Terra wrt {O}
v_SC_i = cell(1,length(dV_tot)); % [km/s] velocità SC prima di entrata SOI Urano wrt {O}
v_i = [ 0 , 1 , 0 ]; % versore velocità SC dopo uscita SOI Terra wrt {O} (ipotizzo tangenziale a traiettoria Terrestre)
v_Urano_f = cell(1,length(dV_tot)); % [km/s] velocità Urano ad encounter
for i = 1:size(dV_tot,2)
    %  velocità di eccesso iperbolica in uscita SOI Terra [km/s]
    [v_inf_T, delta_t, ok, status_msg, info] = deltaV_to_vinf(dV_Terra(1,i),rp_Terra,gm_Terra,SOI_Terra); 
    Earth_exit{i} = struct(...
        'dV_Terra',   dV_Terra(1,i),...  % [km/s] impulso erogato tangenziale a orbita LEO di partenza
        'v_inf_T',    v_inf_T,...        % [km/s] velocità di eccesso iperbolico
        'delta_t',    delta_t,...        % [s] tempo di volo fuga da Terra
        'ok',         ok,...             % [bool] esito fuga
        'status_msg', status_msg,...     % [str] diagnostica
        'info',       info);             % [struct] parametri orbita di fuga (v_c, v_p, v_esc, a, e, F_SOI)
    % vettore velocità SC a uscita SOI Terra wrt {O} [km/s]
    v_SC_i{i} = V_E + v_inf_T * v_i;
    % statovettore SC all'encounter con Urano
    [rx_f,ry_f,rz_f,vx_f,vy_f,vz_f,thetat2,delta_t,ok,status_msg,info] = interplanetary_transfer(r_SC_i(1),r_SC_i(2),r_SC_i(3),v_SC_i{i}(1),v_SC_i{i}(2),v_SC_i{i}(3),gm_Sole,d_Urano);
    r_SC_f{i} = [ rx_f , ry_f , rz_f ];
    v_SC_f{i} = [ vx_f , vy_f , vz_f ];
    X_SC_f{i} = struct( ...
    'r',           [rx_f, ry_f, rz_f], ...  % [km] vettore posizione
    'v',           [vx_f, vy_f, vz_f], ...  % [km/s] vettore velocità
    'theta_t_2',   thetat2, ...             % [rad] anomalia totale = (anomalia vera + argomento al pericentro)
    'delta_t',     delta_t, ...             % [s] tempo di volo della trasferta
    'ok',          ok, ...                  % [bool] esito trasferta
    'status_msg',  status_msg, ...          % [str] diagnostica
    'info',        info );                  % [struct] parametri della trasferta (conic,e,p,a,rp,ra,i,theta1,theta2,branch,rf_chk,delta_t)
    % vettore velocità Urano all'encounter [km/s]
    v_Urano_f{i} = sqrt( gm_Sole / d_Urano ) * [ -sin(thetat2) , cos(thetat2) , 0 ];
    % scalare velocità di eccesso iperbolico Sc in entrata SOI Urano [km/s]
    v_inf(1,i) = norm( v_SC_f{i} - v_Urano_f{i} );
    % determinazione impulso frenata al pericentro orbita iperbolica Urano [km/s]
    [d_V,delta_t,ok,status_msg,info] = Uranus_capture(v_inf(1,i),rp_Urano,gm_Urano,SOI_Urano);
    Uranus_approach{i} = struct(...
        'd_V',        d_V,...            % [km/s] impulso erogato per frenata
        'delta_t',    delta_t,...        % [s] tempo di volo orbita iperbolica di approccio
        'ok',         ok,...             % [bool] esito approccio
        'status_msg', status_msg,...     % [str] diagnostica
        'info',       info );            % [struct] parametri dell’iperbola (a,e,F_SOI, nH, v_p, ecc.)
    if ok
        % determinazione delta_V complessivo
        dV_tot(1,i) = abs(dV_Terra(1,i)) + abs(d_V);
        % determinazione durata missione
        t_durata(1,i) = (Earth_exit{i}.delta_t + X_SC_f{i}.delta_t + Uranus_approach{i}.delta_t) / 86400; % [days] durata della missione
        % determinazione prima finestra utile missione
        t_wait = ( wrapTo2Pi( X_SC_f{i}.theta_t_2 - E_Urano_0 + E_Terra_0 - ( v_Urano / d_Urano - v_Terra / d_Terra ) * Earth_exit{i}.delta_t - ( v_Urano / d_Urano ) * X_SC_f{i}.delta_t ) ) / ( v_Urano / d_Urano ); % [s] tempo che è necessario aspettare per partenza missione
        Et_dep = ET_t0 + t_wait; % partenza in ET
        UTC_departure{i} = cspice_et2utc(Et_dep,'C',2); % partenza in UTC
        Et_arrival = Et_dep + X_SC_f{i}.delta_t + Uranus_approach{i}.delta_t; % arrivo su orbita target in ET
        UTC_arrival{i} = cspice_et2utc(Et_arrival,'C',2); % arrivo su orbita target in UTC
    end  
    SOLUZIONE{i} = struct(...
        'd_V_tot',          dV_tot(1,i),...        % [km/s] impulso erogato totale 
        'data_0',           UTC_departure{i},...   % data di lancio missione (Calendario)
        'data_f',           UTC_arrival{i},...     % data di arrivo orbita target (Calendario)
        'Earth_Exit',       Earth_exit{i},...      % [struct] info su orbita di fuga da Terra
        'Int_transfer',     X_SC_f{i},...          % [struct] info su orbita trasferta interplanetaria
        'Uranus_entrance',  Uranus_approach{i});   % [struct] info su orbita di approccio Urano
end
