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

% definizione range di variabili di progetto (dV) , (r_p) 
% come dV_min viene scelto impulso che comporta trasferta alla Hohmann 
% (se il livello energetico non è già sufficiente dopo flyby) Se invece
% il livello energetico è gia sufficiente dopo flyby dV_min viene scelto
% nullo e dV_max = 5 [km/s]
% come r_p_min viene scelta R_pianeta + h_atmosfera (quota in cui finisce atmosfera)
k1 = 2; % costante moltiplicativa che definisce dV_max = k1 * dV_min
dV_step = 0.1; % [km/s] risoluzione di variabile di progetto dV
k2 = 5; % costante moltiplicativa che definisce r_p_max = k2 * r_p_min
r_p_step = 10; % [km] risoluzione di variabile di progetto r_p


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
SOI_Venere = d_Venere*(gm_Venere/gm_Sole)^5;
SOI_Terra = d_Terra*(gm_Terra/gm_Sole)^5;
SOI_Marte = d_Marte*(gm_Marte/gm_Sole)^5;
SOI_Giove = d_Giove*(gm_Giove/gm_Sole)^5;
SOI_Saturno = d_Saturno*(gm_Saturno/gm_Sole)^5;
SOI_Urano = d_Urano*(gm_Urano/gm_Sole)^5;

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
E_Venere_0 = 2*atan( sqrt(1-elts_Venere_0(2)^2) / (1+elts_Venere_0(2)) * tan(elts_Venere_0(9) / 2) );
E_Terra_0 = 2*atan( sqrt(1-elts_Terra_0(2)^2) / (1+elts_Terra_0(2)) * tan(elts_Terra_0(9) / 2) );
E_Marte_0 = 2*atan( sqrt(1-elts_Marte_0(2)^2) / (1+elts_Marte_0(2)) * tan(elts_Marte_0(9) / 2) );
E_Giove_0 = 2*atan( sqrt(1-elts_Giove_0(2)^2) / (1+elts_Giove_0(2)) * tan(elts_Giove_0(9) / 2) );
E_Saturno_0 = 2*atan( sqrt(1-elts_Saturno_0(2)^2) / (1+elts_Saturno_0(2)) * tan(elts_Saturno_0(9) / 2) );
E_Urano_0 = 2*atan( sqrt(1-elts_Urano_0(2)^2) / (1+elts_Urano_0(2)) * tan(elts_Urano_0(9) / 2) );

%% PRIMA PROPOSTA 
% Terra-->Urano (dV_Terra unica variabile di progetto)

% inizializzazione variabili di progetto 
% dV_Terra [km/s]
rp_Terra = R_Terra + h_LEO; % [km]
Vp_Terra = sqrt( gm_Terra / rp_Terra ); % [km/s]
[dV_Terra] = dV_setter(d_Terra,d_Urano,gm_Sole,gm_Terra,rp_Terra,Vp_Terra,k1,dV_step);

% inizializzazione variabili di costo 
dV_tot = zeros(1,length(dV_Terra));
t_durata = zeros(size(dV_tot));
t_departure = strings(size(dV_tot));

% determinazione velocità di eccesso iperbolica in entrata SOI Urano [km/s]
% inizializzazione
v_inf = zeros(size(dV_tot)); % [km/s] velocità di eccesso iperbolico in entrata SOI Urano
r_SC_f = cell(1,length(dV_tot));
v_SC_f = cell(1,length(dV_tot));
% calcolo
% poichè scelta sdr è arbitraria in questa fase scelgo sdr {O} centrato in Sole
% con x in direzione Terra e z in direzione perpendicolare al piano
% orbitale dei pianeti (ipotesi complanarità)
r_SC_i = d_Terra * [ 1 , 0 , 0 ]; % [km] posizione SC dopo uscita SOI Terra wrt {O}
v_SC_i = cell(1,length(dV_tot)); % [km/s] velocità SC prima di entrata SOI Urano wrt {O}
v_i = [ 0 , 1 , 0 ]; % [km/s] versore velocità SC dopo uscita SOI Terra wrt {O} (ipotizzo tangenziale a traiettoria Terrestre)
v_Urano_f = cell(1,length(dV_tot)); % [km/s] velocità Urano ad encounter
for i = 1:size(dV_tot,2)
    %  velocità di eccesso iperbolica in uscita SOI Terra [km/s]
    v_inf_T = deltaV_to_vinf(dV_Terra(1,i),h_LEO,R_Terra,gm_Terra); 
    % vettore velocità SC a uscita SOI Terra wrt {O} [km/s]
    v_SC_i{i} = v_inf_T * v_i;
    % statovettore SC all'encounter con Urano
    [rx_f,ry_f,rz_f,vx_f,vy_f,vz_f,theta2] = interplanetary_transfer(r_SC_i(1),r_SC_i(2),r_SC_i(3),v_SC_i{i}(1),v_SC_i{i}(2),v_SC_i{i}(3),gm_Sole,d_Urano);
    r_SC_f{i} = [ rx_f , ry_f , rz_f ];
    v_SC_f{i} = [ vx_f , vy_f , vz_f ];
    % vettore velocità Urano all'encounter [km/s]
    v_Urano_f{i} = sqrt( gm_Sole / d_Urano ) * [ -sin(theta2) , cos(theta2) , 0 ];
    % scalare velocità di eccesso iperbolico Sc in entrata SOI Urano [km/s]
    v_inf(1,i) = norm( v_SC_f{i} - v_Urano_f{i} );
end

% determinazione impulso frenata al pericentro orbita iperbolica Urano [km/s]


% determinazione costo proposta 1
% determinazione delta_V complessivo

% determinazione durata missione

% determinazione prima finestra utile missione