%%% MAIN studio preliminare trasferta interplanetaria Terra-Urano
%%% in questo main vengono confrontate rispettivamente X diverse proposte
%%% di missione in X diverse sezioni in cui è suddiviso. Per ogni proposta
%%% vengono analizzate le differenti geometrie di missione che si hanno al
%%% variare dei parametri di progetto (dV_Terra,dV,r_p).
%%% Rispettivamente impulso erogato su orbita LEO di partenza , impulso 
%%% erogato al pericentro dei vari powered-flyby e raggio al pericentro
%%% stesso dei powered-flyby.

% DISCLAIMER: 
% 1) necessaria installazione CSPICE
% 2) Tale studio si basa su metodo delle patched conics
% 3) Tale studio si basa su ipotesi orbite circolari complanari dei pianeti

%% pulizia workspace e command window
clear 
close all
clc
% formattazione calcolo
format long 

%% caricamento Kernels 
cspice_furnsh('kernels\mykernels.furnsh' );

%% Definizione input di missione

% definizione data di riferimento inizio progettazione missione
UTC_t0 = '2025 OCT 30 17:14:37.68'; % t0 in UTC    DISCLAIMER! DA DECIDERE
ET_t0 = cspice_str2et(UTC_t0); % t0 in ET

% definizione orbita LEO di partenza
h_LEO = 1000; % [km] quota orbita LEO circolare di partenza DISCLAIMER da scegliere 

% definizione orbita target intorno Urano
h_target = 50000; % [km] quota orbita circolare target DISCLAIMER da scegliere

%% Definizione range di variabili di progetto (dV_Terra), (dV) , (r_p) 

% dV_Terra
% dV_Terra_min = impulso che comporta trasferta alla Hohmann 
k1 = 2; % costante moltiplicativa che definisce dV_Terra_max = k1 * dV_Terra_min
% dV [km/s] si genera un vettore con step fisso
dV_min = 0;
dV_max = 0; 
dV_step = 0.01; % [km/s] risoluzione di variabile di progetto dV
% r_p [km]  si genera un vettore con step piu denso a quote basse
% r_p_min = R_pianeta + h_atm_pianeta (quota in cui finisce atmosfera densa)
% r_p_step1_pianeta
% r_p_int = R_pianeta + h_int
% r_p_step2_pianeta
% r_p_max = k2_pianeta * R_pianeta



%% estrazione dati kernels su pianeti e lune oggetto di missione


% raggi pianeti e lune [km]
R_Venere = cspice_bodvrd( 'Venus', 'RADII', 3 );
R_Venere = R_Venere(2);  
R_Terra = cspice_bodvrd( 'Earth', 'RADII', 3 );
R_Terra = R_Terra(2);    
R_Marte = cspice_bodvrd( 'Mars', 'RADII', 3 );
R_Marte = R_Marte(2);     
R_Giove = cspice_bodvrd( 'Jupiter', 'RADII', 3 );
R_Giove = R_Giove(2);        
R_Saturno = cspice_bodvrd( 'Saturn', 'RADII', 3 );
R_Saturno = R_Saturno(2);
R_Urano = cspice_bodvrd( 'Uranus', 'RADII', 3 );
R_Urano = R_Urano(2);
R_Titania = cspice_bodvrd( 'Titania', 'RADII', 3 );
R_Titania = R_Titania(2); 
R_Oberon = cspice_bodvrd( 'Oberon', 'RADII', 3 );
R_Oberon = R_Oberon(2);  

% parametri gravitazionali pianeti e lune [km^3/s^2]
gm_Sole = cspice_bodvrd('Sun','GM',1);
gm_Venere = cspice_bodvrd('Venus','GM',1);
gm_Terra = cspice_bodvrd('Earth','GM',1);
gm_Marte = cspice_bodvrd('Mars','GM',1);
gm_Giove = cspice_bodvrd('Jupiter','GM',1);
gm_Saturno = cspice_bodvrd('Saturn','GM',1);
gm_Urano = cspice_bodvrd('Uranus','GM',1);
gm_Titania = cspice_bodvrd('Titania','GM',1);
gm_Oberon = cspice_bodvrd('Oberon','GM',1);

% raggi orbite circolari pianeti [km]
d_Venere =   108000000;
d_Terra =    149598000;
d_Marte =    228000000;
d_Giove =    778412020;
d_Saturno =  1427000000;
d_Urano =    2870660000;

% raggi SOI dei pianeti [km]
SOI_Venere  = d_Venere  * (gm_Venere/gm_Sole)^(2/5);
SOI_Terra   = d_Terra   * (gm_Terra /gm_Sole)^(2/5);
SOI_Marte   = d_Marte   * (gm_Marte /gm_Sole)^(2/5);
SOI_Giove   = d_Giove   * (gm_Giove /gm_Sole)^(2/5);
SOI_Saturno = d_Saturno * (gm_Saturno/gm_Sole)^(2/5);
SOI_Urano   = d_Urano   * (gm_Urano /gm_Sole)^(2/5);

% velocità pianeti [km/s]
v_Venere = sqrt(gm_Sole/d_Venere);
v_Terra = sqrt(gm_Sole/d_Terra);
v_Marte = sqrt(gm_Sole/d_Marte);
v_Giove = sqrt(gm_Sole/d_Giove);
v_Saturno = sqrt(gm_Sole/d_Saturno);
v_Urano = sqrt(gm_Sole/d_Urano);

%% determinazione condizioni iniziali

% estrazione stato vettore dei pianeti al tempo t0 in J2000 centrato in Sole [km , km/s]
[x_Venere_0,ltimeV]= cspice_spkezr('Venus',ET_t0,'J2000','NONE','Sun');
[x_Terra_0,ltimeT]= cspice_spkezr('Earth',ET_t0,'J2000','NONE','Sun');
[x_Marte_0,ltimeM]= cspice_spkezr('4',ET_t0,'J2000','NONE','Sun');
[x_Giove_0,ltimeG]= cspice_spkezr('5',ET_t0,'J2000','NONE','Sun');
[x_Saturno_0,ltimeS]= cspice_spkezr('6',ET_t0,'J2000','NONE','Sun');
[x_Urano_0,ltimeU]= cspice_spkezr('7',ET_t0,'J2000','NONE','Sun');
% NOTA BENE per i pianeti più esterni rispetto alla Terra si estrae lo
% stato vettore del baricentro del sistema pianeta + satelliti poichè non
% sono disponibili kernel che si riferiscano al singolo pianeta

% determinazione anomalia eccentrica pianeti al tempo t0 rispetto a J2000 centrato nel Sole [rad]
%1) estrazione parametri orbitali 
elts_Venere_0 = cspice_oscltx(x_Venere_0,ET_t0,gm_Sole);
elts_Terra_0 = cspice_oscltx(x_Terra_0,ET_t0,gm_Sole);
elts_Marte_0 = cspice_oscltx(x_Marte_0,ET_t0,gm_Sole);
elts_Giove_0 = cspice_oscltx(x_Giove_0,ET_t0,gm_Sole);
elts_Saturno_0 = cspice_oscltx(x_Saturno_0,ET_t0,gm_Sole);
elts_Urano_0 = cspice_oscltx(x_Urano_0,ET_t0,gm_Sole);
%2) trasformazione true anomaly in eccentric anomaly [rad]
E_Venere_0 = wrapTo2Pi( 2*atan( sqrt(1-elts_Venere_0(2)^2) / (1+elts_Venere_0(2)) * tan(elts_Venere_0(9) / 2) ) );
E_Terra_0 = wrapTo2Pi( 2*atan( sqrt(1-elts_Terra_0(2)^2) / (1+elts_Terra_0(2)) * tan(elts_Terra_0(9) / 2) ) );
E_Marte_0 = wrapTo2Pi( 2*atan( sqrt(1-elts_Marte_0(2)^2) / (1+elts_Marte_0(2)) * tan(elts_Marte_0(9) / 2) ) );
E_Giove_0 = wrapTo2Pi( 2*atan( sqrt(1-elts_Giove_0(2)^2) / (1+elts_Giove_0(2)) * tan(elts_Giove_0(9) / 2) ) );
E_Saturno_0 = wrapTo2Pi( 2*atan( sqrt(1-elts_Saturno_0(2)^2) / (1+elts_Saturno_0(2)) * tan(elts_Saturno_0(9) / 2) ) );
E_Urano_0 = wrapTo2Pi( 2*atan( sqrt(1-elts_Urano_0(2)^2) / (1+elts_Urano_0(2)) * tan(elts_Urano_0(9) / 2) ) );

%% determinazione raggi di orbite circolari di partenza e target [km]

rp_Terra = R_Terra + h_LEO; 
rp_Urano = R_Urano + h_target; 

%% PRIMA PROPOSTA 
% Terra-->Urano (dV_Terra) variabile di progetto

% velocità SC in orbita LEO [km/s]
Vp_Terra = sqrt( gm_Terra / rp_Terra ); 
% inizializzazione variabili di progetto 
% dV_Terra [km/s]
[dV_Terra] = dV_setter(d_Terra,d_Urano,gm_Sole,gm_Terra,rp_Terra,Vp_Terra,k1,dV_max,dV_step);

% soluzioni prima proposta
[SOLUZIONE_T_U] = TERRA_URANO(d_Terra,d_Urano,gm_Terra,gm_Urano,gm_Sole,SOI_Terra,SOI_Urano,v_Terra,v_Urano,rp_Terra,rp_Urano,E_Terra_0,E_Urano_0,ET_t0,dV_Terra);

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

%% SECONDA PROPOSTA
% Terra-->Marte-->Giove-->Urano
% (dV_Terra,dV_Marte,rp_Marte,dV_Giove,rp_Giove) variabili di progetto

% inizializzazione variabili di progetto
% dV_Terra [km/s]
[dV_Terra] = dV_setter(d_Terra,d_Marte,gm_Sole,gm_Terra,rp_Terra,Vp_Terra,k1,dV_max,dV_step);
% dV [km/s]
dV = dV_min:dV_step:dV_max;
% r_p [km]
h_atm_Marte = 200; % quota atmosfera [km] DA CAMBIARE
h_int_Marte = 2000; % [km] quota intermedia
k2_Marte = 4;
r_p_step_Marte1 = 200;
r_p_step_Marte2 = 1000;
h_atm_Giove = 4000; % quota atmosfera Giove [km] DA CAMBIARE
h_int_Giove = R_Giove;
k2_Giove = 20;
r_p_step_Giove1 = 0.2 * R_Giove;
r_p_step_Giove2 = 0.5 * R_Giove;
r_p_Marte = unique([R_Marte + h_atm_Marte : r_p_step_Marte1 : R_Marte + h_int_Marte , R_Marte + h_int_Marte : r_p_step_Marte2 : R_Marte * k2_Marte]) ;
r_p_Giove = unique([R_Giove + h_atm_Giove : r_p_step_Giove1 : R_Giove + h_int_Giove , R_Giove + h_int_Giove : r_p_step_Giove2 : R_Giove * k2_Giove]) ;


[SOLUZIONE_T_M_G_U] = TERRA_MARTE_GIOVE_URANO(d_Terra,d_Marte,d_Giove,d_Urano,gm_Terra,gm_Marte,gm_Giove,gm_Urano,gm_Sole,SOI_Terra,SOI_Marte,SOI_Giove,SOI_Urano,v_Terra,v_Marte,v_Giove,v_Urano,rp_Terra,rp_Urano,E_Terra_0,E_Marte_0,E_Giove_0,E_Urano_0,ET_t0,dV_Terra,dV,r_p_Marte,r_p_Giove);
[best_soluzione, idx] = pick_best_solution(SOLUZIONE_T_M_G_U);
disp(best_soluzione)    % la struct della cella migliore
disp(idx)               % indici cella migliore
