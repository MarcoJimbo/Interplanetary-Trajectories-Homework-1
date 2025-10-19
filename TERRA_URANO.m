function[SOLUZIONE] = TERRA_URANO(d_Terra,d_Urano,gm_Terra,gm_Urano,gm_Sole,SOI_Terra,SOI_Urano,v_Terra,v_Urano,rp_Terra,rp_Urano,E_Terra_0,E_Urano_0,ET_t0,dV_Terra)
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
%%% SOLUZIONE:                                     CELL   [1xn] of struct of possible solutions performances
%%%          .d_V_tot:                             SCALAR [1x1] = impulso totale missione [km/s]
%%%          .data_0:                              STRING       = data di lancio (calendario)
%%%          .data_f:                              STRING       = data di arrivo (calendario)
%%%          .Earth_Exit:                          STRUCT [1x1] = info orbita di fuga Terra
%%%                     .dV_terra:                 SCALAR [1x1] = impulso erogato su LEO [km/s]
%%%                     .v_inf_T:                  SCALAR [1x1] = velocità di eccesso iperbolico [km/s]
%%%                     .delta_t:                  SCALAR [1x1] = tempo di volo [s]
%%%                     .ok:                       BOOL         = esito fuga
%%%                     .status_msg:               STRING       = diagnostica
%%%                     .info:                     STRUCT [1x1] = info orbita
%%%                          .branch:              STRING       = conic
%%%                          .r_p:                 SCALAR [1x1] = raggio al pericentro [km]
%%%                          .v_c:                 SCALAR [1x1] = velocità LEO [km/s]
%%%                          .v_p:                 SCALAR [1x1] = velocità dopo impulso [km/s]
%%%                          .v_esc:               SCALAR [1x1] = velocità di fuga parabolica [km/s]
%%%                          .a:                   SCALAR [1x1] = semi asse maggiore [km]
%%%                          .e:                   SCALAR [1x1] = eccentricità
%%%                          .F_SOI:               SCALAR [1x1] = anomalia iperbolica al raggio della SOI [rad]
%%%                          .nH:                  SCALAR [1x1] = anomalia media iperbolica [rad/s]
%%%                          .r_SOI:               SCALAR [1x1] = raggio SOI [km]
%%%          .Int_Transfer
%%%                       .r:                      DOUBLE [1x3] = posizione SC ad encounter Urano wrt s.d.r. {O} [km]
%%%                       .v:                      DOUBLE [1x3] = velocità SC ad encounter Urano wrt s.d.r {O} [km/s]
%%%                       .theta_t_2:              SCALAR [1x1] = anomalia totale ad encounter wrt s.d.r {O} [rad]
%%%                       .delta_t:                SCALAR [1x1] = tempo di volo [s]
%%%                       .ok:                     BOOL         = esito
%%%                       .status_msg:             STRING       = diagnostica
%%%                       .info                    STRUCT [1x1] = info orbita interplanetaria
%%%                            .conic:             STRING       = conic type
%%%                            .e:                 SCALAR [1x1] = eccentricità
%%%                            .p:                 SCALAR [1x1] = semilato retto [km]
%%%                            .a:                 SCALAR [1x1] = semiasse maggiore [km]
%%%                            .rp:                SCALAR [1x1] = raggio al pericentro [km]
%%%                            .ra:                SCALAR [1x1] = raggio all'apocentro [km]
%%%                            .i:                 SCALAR [1x1] = inclinazione [rad]
%%%                            .theta1:            SCALAR [1x1] = anomalia vera iniziale [rad]
%%%                            .theta2:            SCALAR [1x1] = anomalia vera encounter [rad]
%%%                            .branch:            STRING       = direzione trasferta
%%%                            .rf_chk:            SCALAR [1x1] = check rf raggiunto [km]
%%%                            .delta_t:           SCALAR [1x1] = tempo di volo [s]
%%%          .Uranus_Entrance
%%%                          .d_V:                 SCALAR [1x1] = impulso frenata su Urano [km/s]
%%%                          .delta_t:             SCALAR [1x1] = tempo di volo [s]
%%%                          .ok:                  BOOL         = esito
%%%                          .status_msg:          STRING       = diagnostica
%%%                          .info:                STRUCT [1x1] = info orbita
%%%                               .a:              SCALAR [1x1] = semiasse maggiore [km]
%%%                               .e:              SCALAR [1x1] = eccentricità
%%%                               .F_SOI:          SCALAR [1x1] = anomalia iperbolica al raggio della SOI [rad]
%%%                               .nH:             SCALAR [1x1] = anomalia media iperbolica [rad/s]
%%%                               .v_p:            SCALAR [1x1] = velocità al pericentro iperbole [km/s]
%%%                               .v_c:            SCALAR [1x1] = velocità circolare orbita target [km/s]
%%%                               .r_SOI:          SCALAR [1x1] = raggio SOI [km]


% inizializzazione variabili di costo 
dV_tot = NaN(1,length(dV_Terra));
t_durata = zeros(size(dV_tot));
UTC_departure = strings(size(dV_tot));
UTC_arrival = strings(size(dV_tot));


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
        w_U = v_Urano / d_Urano; % [rad/s]
        w_T = v_Terra / d_Terra; % [rad/s]
        % inizializzazione ciclo while 
        k = 0;
        t_wait = -1;
        while t_wait < 0
        t_wait = ( E_Terra_0 - E_Urano_0 - w_U * Earth_exit{i}.delta_t - w_U * X_SC_f{i}.delta_t + w_T * Earth_exit{i}.delta_t + thetat2 ) / ( w_U - w_T ) + ( 2 * k * pi ) / abs( w_U - w_T ); % [s] tempo che è necessario aspettare per partenza missione
        k = k + 1;
        end
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
