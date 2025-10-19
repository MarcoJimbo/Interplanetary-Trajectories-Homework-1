function[SOLUZIONE] = TERRA_MARTE_GIOVE_URANO(d_Terra,d_Marte,d_Giove,d_Urano,gm_Terra,gm_Marte,gm_Giove,gm_Urano,gm_Sole,SOI_Terra,SOI_Marte,SOI_Giove,SOI_Urano,v_Terra,v_Marte,v_Giove,v_Urano,rp_Terra,rp_Urano,E_Terra_0,E_Marte_0,E_Giove_0,E_Urano_0,ET_t0,dV_Terra,dV,r_p_Marte,r_p_Giove)
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
dV_tot = NaN(length(dV_Terra),length(dV),length(r_p_Marte),length(dV),length(r_p_Giove));
t_durata = zeros(size(dV_tot));
UTC_departure = strings(size(dV_tot));
UTC_arrival = strings(size(dV_tot));

% inizializzazione
SOLUZIONE = cell(size(dV_tot)); % cella contenente tutte le caratteristiche della missione
X_Mars_encounter = cell(1,size(dV_tot,1));
X_Mars_exit = cell(size(dV_tot,1),size(dV_tot,2),size(dV_tot,3));
X_Jupiter_encounter = cell(size(dV_tot,1),size(dV_tot,2),size(dV_tot,3));
X_Jupiter_exit = cell(size(dV_tot));
X_Uranus_encounter = cell(size(dV_tot)); % statovettore dello Spacecraft all'encounter + informazioni su traiettoria interplanetaria
v_inf_U = zeros(size(dV_tot)); % [km/s] velocità di eccesso iperbolico in entrata SOI Urano
Earth_exit = cell(1,size(dV_tot,1)); % informazioni su traiettoria fuga da Terra
Mars_flyby = cell(size(dV_tot,1),size(dV_tot,2),size(dV_tot,3)); % informazioni su traiettoria flyby Marte
Jupiter_flyby = cell(size(dV_tot)); % informazioni su traiettoria flyby Giove
Uranus_approach = cell(size(dV_tot)); % informazioni su traiettoria iperbolica in SOI Urano 
V_E = [ 0 , v_Terra , 0];
% calcolo
% poichè scelta sdr è arbitraria in questa fase scelgo sdr inerziale {O} centrato in Sole
% con x in direzione Terra (al momento SC esce da SOI Terra) e z in direzione perpendicolare al piano
% orbitale dei pianeti (ipotesi complanarità)
r_SC_i = d_Terra * [ 1 , 0 , 0 ]; % [km] posizione SC dopo uscita SOI Terra wrt {O}
v_SC_i = cell(1,size(dV_tot,1)); % [km/s] velocità SC dopo uscita SOI Terra wrt {O}
v_i = [ 0 , 1 , 0 ]; % versore velocità SC dopo uscita SOI Terra wrt {O} (ipotizzo tangenziale a traiettoria Terrestre)
v_Marte_f = cell(1,size(dV_tot,1)); % [km/s] velocità Marte ad encounter
v_Giove_f = cell(size(dV_tot,1),size(dV_tot,2),size(dV_tot,3)); % [km/s] velocità Giove ad encounter
v_Urano_f = cell(size(dV_tot)); % [km/s] velocità Urano ad encounter

for i = 1:size(dV_tot,1) % variazione di dV_Terra
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

    % statovettore SC all'encounter con Marte
    [rx_f,ry_f,rz_f,vx_f,vy_f,vz_f,thetat2,delta_t,ok,status_msg,info] = interplanetary_transfer(r_SC_i(1),r_SC_i(2),r_SC_i(3),v_SC_i{i}(1),v_SC_i{i}(2),v_SC_i{i}(3),gm_Sole,d_Marte);
    X_Mars_encounter{i} = struct( ...
    'r',           [rx_f, ry_f, rz_f], ...  % [km] vettore posizione
    'v',           [vx_f, vy_f, vz_f], ...  % [km/s] vettore velocità
    'theta_t_2',   thetat2, ...             % [rad] anomalia totale = (anomalia vera + argomento al pericentro)
    'delta_t',     delta_t, ...             % [s] tempo di volo della trasferta
    'ok',          ok, ...                  % [bool] esito trasferta
    'status_msg',  status_msg, ...          % [str] diagnostica
    'info',        info );                  % [struct] parametri della trasferta (conic,e,p,a,rp,ra,i,theta1,theta2,branch,rf_chk,delta_t)
    % vettore velocità Marte all'encounter [km/s]
    v_Marte_f{i} = sqrt( gm_Sole / d_Marte ) * [ -sin(thetat2) , cos(thetat2) , 0 ];
    % scalare velocità di eccesso iperbolico Sc in entrata SOI Marte [km/s]
    v_inf_M = X_Mars_encounter{i}.v - v_Marte_f{i} ;
    norm_v_inf_M = norm(v_inf_M);

    % statovettore SC all'uscita da SOI Marte
    for j=1:size(dV_tot,2)
        for l=1:size(dV_tot,3)
            [v_inf2, def_angle, t_in, t_out, delta_t, ok, status_msg, info] = general_flyby(norm_v_inf_M, r_p_Marte(1,l), gm_Marte, dV(1,j), SOI_Marte);
            r_SC = X_Mars_encounter{i}.r;
            v_SC = ([cos(def_angle),sin(def_angle),0 ; -sin(def_angle),cos(def_angle),0 ; 0,0,1]' * (v_inf_M*v_inf2/norm_v_inf_M)')' + v_Marte_f{i};
            X_Mars_exit{i,j,l} = struct( ...
            'rp',          r_p_Marte(1,l),...      % [km] raggio al pericentro
            'dV',          dV(1,j),...             % [km/s] impulso erogato
            'r',           [r_SC(1) , r_SC(2), r_SC(3)], ...  % [km] vettore posizione
            'v',           [v_SC(1) , v_SC(2), v_SC(3)], ...  % [km/s] vettore velocità
            'v_inf_1',     norm_v_inf_M, ...        % [km/s] velocità di eccesso iperbolico in entrata
            't_in',        t_in,...                 % [s] tempo di volo pre manovra
            't_out',       t_out,...                % [s] tempo di volo post manovra
            'delta_t',     delta_t, ...             % [s] tempo di volo della trasferta
            'ok',          ok, ...                  % [bool] esito trasferta
            'status_msg',  status_msg, ...          % [str] diagnostica
            'info',        info );                  % [struct] parametri della trasferta (conic,e,p,a,rp,ra,i,theta1,theta2,branch,rf_chk,delta_t)
            
            % statovettore SC all'encounter con Giove
            [rx_f,ry_f,rz_f,vx_f,vy_f,vz_f,thetat2,delta_t,ok,status_msg,info] = interplanetary_transfer(r_SC(1),r_SC(2),r_SC(3),v_SC(1),v_SC(2),v_SC(3),gm_Sole,d_Giove);
            r_SC = [ rx_f , ry_f , rz_f ];
            v_SC = [ vx_f , vy_f , vz_f ];
            X_Jupiter_encounter{i,j,l} = struct( ...
            'r',           [rx_f, ry_f, rz_f], ...  % [km] vettore posizione
            'v',           [vx_f, vy_f, vz_f], ...  % [km/s] vettore velocità
            'theta_t_2',   thetat2, ...             % [rad] anomalia totale = (anomalia vera + argomento al pericentro)
            'delta_t',     delta_t, ...             % [s] tempo di volo della trasferta
            'ok',          ok, ...                  % [bool] esito trasferta
            'status_msg',  status_msg, ...          % [str] diagnostica
            'info',        info );                  % [struct] parametri della trasferta (conic,e,p,a,rp,ra,i,theta1,theta2,branch,rf_chk,delta_t)
            % vettore velocità Giove all'encounter [km/s]
            v_Giove_f{i,j,l} = sqrt( gm_Sole / d_Giove ) * [ -sin(thetat2) , cos(thetat2) , 0 ];
            % scalare velocità di eccesso iperbolico Sc in entrata SOI Urano [km/s]
            v_inf_G =  v_SC - v_Giove_f{i,j,l} ;
            norm_v_inf_G = norm(v_inf_G);

            % statovettore SC all'uscita da SOI Giove
            for k=1:size(dV_tot,4)
                for m=1:size(dV_tot,5)

                    [v_inf2, def_angle, t_in, t_out, delta_t, ok, status_msg, info] = general_flyby(norm_v_inf_G, r_p_Giove(1,m), gm_Giove , dV(1,k), SOI_Giove);
                    r_SC = X_Jupiter_encounter{i}.r;
                    v_SC = ([cos(def_angle),sin(def_angle),0 ; -sin(def_angle),cos(def_angle),0 ; 0,0,1]' * (v_inf_G*v_inf2/norm_v_inf_G)')' + v_Giove_f{i,j,l};
                    X_Jupiter_exit{i,j,l,k,m} = struct( ...
                    'rp',          r_p_Giove(1,m),...      % [km] raggio al pericentro
                    'dV',          dV(1,k),...             % [km/s] impulso erogato
                    'r',           [r_SC(1) , r_SC(2), r_SC(3)], ...  % [km] vettore posizione
                    'v',           [v_SC(1) , v_SC(2), v_SC(3)], ...  % [km/s] vettore velocità
                    'v_inf_1',     norm_v_inf_G, ...        % [km/s] velocità di eccesso iperbolico in entrata
                    't_in',        t_in,...                 % [s] tempo di volo pre manovra
                    't_out',       t_out,...                % [s] tempo di volo post manovra
                    'delta_t',     delta_t, ...             % [s] tempo di volo della trasferta
                    'ok',          ok, ...                  % [bool] esito trasferta
                    'status_msg',  status_msg, ...          % [str] diagnostica
                    'info',        info );                  % [struct] parametri della trasferta (conic,e,p,a,rp,ra,i,theta1,theta2,branch,rf_chk,delta_t)
                    
                    % statovettore SC all'encounter con Urano
                    [rx_f,ry_f,rz_f,vx_f,vy_f,vz_f,thetat2,delta_t,ok,status_msg,info] = interplanetary_transfer(r_SC(1),r_SC(2),r_SC(3),v_SC(1),v_SC(2),v_SC(3),gm_Sole,d_Urano);
                    r_SC = [ rx_f , ry_f , rz_f ];
                    v_SC = [ vx_f , vy_f , vz_f ];
                    X_Uranus_encounter{i,j,l,k,m} = struct( ...
                    'r',           [rx_f, ry_f, rz_f], ...  % [km] vettore posizione
                    'v',           [vx_f, vy_f, vz_f], ...  % [km/s] vettore velocità
                    'theta_t_2',   thetat2, ...             % [rad] anomalia totale = (anomalia vera + argomento al pericentro)
                    'delta_t',     delta_t, ...             % [s] tempo di volo della trasferta
                    'ok',          ok, ...                  % [bool] esito trasferta
                    'status_msg',  status_msg, ...          % [str] diagnostica
                    'info',        info );                  % [struct] parametri della trasferta (conic,e,p,a,rp,ra,i,theta1,theta2,branch,rf_chk,delta_t)
                    % vettore velocità Urano all'encounter [km/s]
                    v_Urano_f{i,j,l,k,m} = sqrt( gm_Sole / d_Urano ) * [ -sin(thetat2) , cos(thetat2) , 0 ];
                    % scalare velocità di eccesso iperbolico Sc in entrata SOI Urano [km/s]
                    v_inf_U =  v_SC - v_Urano_f{i,j,l,k,m} ;
                    norm_v_inf_U = norm(v_inf_U);

                    % determinazione impulso frenata al pericentro orbita iperbolica Urano [km/s]
                    [d_V,delta_t,ok,status_msg,info] = Uranus_capture(norm_v_inf_U,rp_Urano,gm_Urano,SOI_Urano);
                    Uranus_approach{i,j,l,k,m} = struct(...
                    'd_V',        d_V,...            % [km/s] impulso erogato per frenata
                    'delta_t',    delta_t,...        % [s] tempo di volo orbita iperbolica di approccio
                    'ok',         ok,...             % [bool] esito approccio
                    'status_msg', status_msg,...     % [str] diagnostica
                    'info',       info );            % [struct] parametri dell’iperbola (a,e,F_SOI, nH, v_p, ecc.)
                    if ok
                    % determinazione delta_V complessivo
                    dV_tot(i,j,l,k,m) = abs(dV_Terra(1,i)) + abs(dV(1,j)) + abs(dV(1,k)) + abs(Uranus_approach{i,j,l,k,m}.d_V);
                    % determinazione durata missione
                    t_durata(i,j,l,k,m) = (Earth_exit{i}.delta_t + X_Mars_encounter{i}.delta_t + X_Mars_exit{i,j,l}.delta_t + X_Jupiter_encounter{i,j,l}.delta_t + X_Jupiter_exit{i,j,l,k,m}.delta_t + X_Uranus_encounter{i,j,l,k,m}.delta_t + Uranus_approach{i,j,l,k,m}.delta_t) / 86400; % [days] durata della missione
                    % determinazione prima finestra utile missione
                    % velocità angolari pianeti [rad/s]
                    w_T = v_Terra / d_Terra; 
                    w_M = v_Marte / d_Marte;
                    w_G = v_Giove / d_Giove;
                    w_U = v_Urano / d_Urano;
                    % periodi sinodici % [s]
                    T_TM = ( 2*pi ) / abs( w_T - w_M );
                    T_TG = ( 2*pi ) / abs( w_T - w_G );
                    T_TU = ( 2*pi ) / abs( w_T - w_U );
                    % tempo di attesa prima di lancio [s]
                    t_wait1 = ( X_Mars_encounter{i}.theta_t_2 - E_Marte_0 + E_Terra_0 - w_M * ( Earth_exit{i}.delta_t + X_Mars_encounter{i}.delta_t ) + w_T * Earth_exit{i}.delta_t ) / ( w_M - w_T ); % [s] 
                    t_wait2 = ( X_Jupiter_encounter{i,j,l}.theta_t_2 - E_Giove_0 + E_Terra_0 - w_G * ( Earth_exit{i}.delta_t + X_Mars_encounter{i}.delta_t + X_Mars_exit{i,j,l}.delta_t + X_Jupiter_encounter{i,j,l}.delta_t ) + w_T * Earth_exit{i}.delta_t ) / ( w_G - w_T );
                    t_wait3 = ( X_Uranus_encounter{i,j,l}.theta_t_2 - E_Urano_0 + E_Terra_0 - w_U * ( Earth_exit{i}.delta_t + X_Mars_encounter{i}.delta_t + X_Mars_exit{i,j,l}.delta_t + X_Jupiter_encounter{i,j,l}.delta_t + X_Jupiter_exit{i,j,l,k,m}.delta_t + X_Uranus_encounter{i,j,l,k,m}.delta_t ) + w_T * Earth_exit{i}.delta_t ) / ( w_U - w_T );
                    epsTol_days = 15; % [s]
                    epsTol = epsTol_days*86400;
                    dt_user = epsTol/5;
                    [t_wait, T_repeat, det] = common_wait_solver( [t_wait1;t_wait2;t_wait3] , [T_TM;T_TG;T_TU] , epsTol ,dt_user);
                    % if det.success
                    %     fprintf('t_min = %.3f giorni\n', t_wait/86400);
                    %     fprintf('T_repeat ≈ %.3f anni\n', T_repeat/86400/365.25);
                    %     disp('Residui (s):'), disp(det.residuals)
                    %     disp(det.method), fprintf('Δt usato = %.2f h\n', det.dt/3600);
                    % else
                    %     fprintf('Nessuna tripla soluzione entro l’orizzonte (%.1e anni).\n', det.extra.Kmax*det.dt*det.m12/86400/365.25);
                    %     fprintf('Miglior scarto trovato: %.1f s (≈ %.2f giorni) a t ≈ %.3f anni\n', ...
                    %         det.extra.best_residual_sec, det.extra.best_residual_sec/86400, det.extra.best_t_sec/86400/365.25);
                    % end
                    Et_dep = ET_t0 + t_wait; % partenza in ET
                    UTC_departure{i,j,l,k,m} = cspice_et2utc(Et_dep,'C',2); % partenza in UTC
                    Et_arrival = Et_dep + X_Uranus_encounter{i,j,l,k,m}.delta_t + Uranus_approach{i,j,l,k,m}.delta_t; % arrivo su orbita target in ET
                    UTC_arrival{i,j,l,k,m} = cspice_et2utc(Et_arrival,'C',2); % arrivo su orbita target in UTC
                    end  
                    SOLUZIONE{i,j,l,k,m} = struct(...
                    'd_V_tot',           dV_tot(i,j,l,k,m),...            % [km/s] impulso erogato totale 
                    'data_0',            UTC_departure{i,j,l,k,m},...     % data di lancio missione (Calendario)
                    'data_f',            UTC_arrival{i,j,l,k,m},...       % data di arrivo orbita target (Calendario)
                    'Earth_Exit',        Earth_exit{i},...                % [struct] info su orbita di fuga da Terra
                    'Earth_Mars_tr',     X_Mars_encounter{i},...          % [struct] info su orbita trasferta interplanetaria
                    'Mars_flyby',        X_Mars_exit{i,j,l},...
                    'Mars_Jupiter_tr',   X_Jupiter_encounter{i,j,l},...   % [struct] info su orbita trasferta interplanetaria
                    'Jupiter_flyby',     X_Jupiter_exit{i,j,l,k,m},...
                    'Jupiter_Uranus_tr', X_Uranus_encounter{i,j,l,k,m},...% [struct] info su orbita trasferta interplanetaria
                    'Uranus_entrance',   Uranus_approach{i,j,l,k,m});     % [struct] info su orbita di approccio Urano
                end
            end
        end
    end

end
