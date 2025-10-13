function[rx_f,ry_f,rz_f,vx_f,vy_f,vz_f,theta_t_2] = interplanetary_transfer(rx_i,ry_i,rz_i,vx_i,vy_i,vz_i,gm_Sole,rf)
%%% funzione che dato in input lo statovettore dello Spacecraft e il parametro 
%%% gravitazionale del Sole ne  determina lo statovettore quando si 
%%% trova a distanza rf da esso. DISCLAIMER si fornisce solo il primo dei
%%% due statovettori soluzione che soddisfano r = rf!

%%% INPUT:
%%% 1) rx_i = [km] componente x della posizione iniziale dello SC in generico s.d.r. centrato nel Sole
%%% 1) ry_i = [km] componente y della posizione iniziale dello SC in generico s.d.r. centrato nel Sole
%%% 1) rz_i = [km] componente z della posizione iniziale dello SC in generico s.d.r. centrato nel Sole
%%% 1) vx_i = [km/s] componente x della velocità iniziale dello SC in generico s.d.r. centrato nel Sole
%%% 1) vy_i = [km/s] componente y della velocità iniziale dello SC in generico s.d.r. centrato nel Sole
%%% 1) vz_i = [km/s] componente z della velocità iniziale dello SC in generico s.d.r. centrato nel Sole
%%% 1) gm_Sole = [km^3/s^2] parametro gravitazionale Sole
%%% 1) rf = [km] distanza dal Sole a cui si vuole determinare statovettore

%%% OUTPUT:
%%% 1) rx_f = [km] componente x della posizione finale dello SC in generico s.d.r. centrato nel Sole
%%% 2) ry_f = [km] componente y della posizione finale dello SC in generico s.d.r. centrato nel Sole
%%% 3) rz_f = [km] componente z della posizione finale dello SC in generico s.d.r. centrato nel Sole
%%% 4) vx_f = [km/s] componente x della velocità finale dello SC in generico s.d.r. centrato nel Sole
%%% 5) vy_f = [km/s] componente y della velocità finale dello SC in generico s.d.r. centrato nel Sole
%%% 1) vz_f = [km/s] componente z della velocità finale dello SC in generico s.d.r. centrato nel Sole
%%% 7) theta_t_2 = [rad] anomalia vera + argomento al pericentro di encounter



% vettore posizione uscita SOI pianeta di partenza wrt J2000 centrato in Sole [km]
r_i = [rx_i ; ry_i ; rz_i]; 

% vettore velocità eliocentrica uscita SOI pianeta di partenza wrt J2000 centrato in Sole [km/s] 
v_i = [vx_i ; vy_i ; vz_i];

% calcolo parametri orbitali 
r = norm(r_i);
v = norm(v_i);
energy = -r / gm_Sole + v^2 / 2; % [km^2/s^2] energia orbita eliocentrica
a = - gm_Sole / (2*energy); % [km] semiasse maggiore
h_vec = cross(r_i,v_i); % [km^2/s]
h = norm(h_vec);
p = h^2 / gm_Sole; % [km]
e = sqrt( 1 - p / a );
vr = dot(v_i,r_i) / r; % [km/s]
So1 = vr / e * sqrt( a*(1-e^2) / gm_Sole);
Co1 = ( p/r - 1 ) / e;
theta1 = 2 * atan( So1 / (1-Co1) ); % [rad] anomalia vera 
vt = sqrt( gm_Sole / p ) * e * sin(theta1);
r_versore = r_i / r;
h_versore = h_vec / h;
theta_versore = cross(h_versore,r_versore) / norm( cross(h_versore,r_versore) );
R_LVLH_ECI_t0 = [ r_versore , theta_versore , h_versore];
i = acos(h_versore(3)); % [rad] inclinazione orbita eliocentica
% tolleranza su errore numerico
if i < 1e-2
    i = 0;
end
if i == 0 % casistica traiettoria planare
    Cot1 = R_LVLH_ECI_t0(1,1);
    Sot1 = R_LVLH_ECI_t0(2,1);
    theta_t_1 = 2 * atan( Sot1 / ( 1 + Cot1 ) ); % [rad]
    w = theta_t_1 - theta1; % [rad] 
    % OMEGA non definito
else % casistica traiettoria 3D
    Sot1 = R_LVLH_ECI_t0(3,1) / sin(i);
    Cot1 = R_LVLH_ECI_t0(3,2) / sin(i);
    theta_t_1 = 2 * atan( Sot1 / ( 1 + Cot1 ) ); % [rad]
    w = theta_t_1 - theta1; % [rad]
    SOM1 = R_LVLH_ECI_t0(1,3) / sin(i);
    COM1 = -R_LVLH_ECI_t0(2,3) / sin(i);
    OM_1 = 2 * atan( SOM1 / (1+COM1) ); % [rad]
end

% determinazione anomalia vera al primo encounter [rad]
if v > sqrt( gm_Sole / r ) % trasferta verso pianeta esterno
    Cot2 = ( p/rf - 1 ) / e;
    if Cot2 < -1 || Cot2 > 1
        error('La traiettoria eliocentrica non raggiunge la distanza dal sole rf specificata in input')
    end
    theta21 = acos( ( p/rf - 1 ) / e ); 
    theta22 = 2*pi - acos( ( p/rf - 1 ) / e );
    vr21 = sqrt( gm_Sole / p ) * e * sin(theta21);
    vr22 = sqrt( gm_Sole / p ) * e * sin(theta22);
    if vr21 > 0 && vr22 < 0
        theta2 = theta21;
    elseif  vr21 < 0 && vr22 > 0
        theta2 = theta22;
    else 
        error('errore nella risoluzione del calcolo dell anomalia vera al primo encounter')
    end
elseif v < sqrt( gm_Sole / r ) % trasferta verso pianeta interno
    Cot2 = ( p/rf - 1 ) / e;
    if Cot2 < -1 || Cot2 > 1
        error('La traiettoria eliocentrica non raggiunge la distanza dal sole rf specificata in input')
    end
    theta21 = acos( ( p/rf - 1 ) / e ); 
    theta22 = 2*pi - acos( ( p/rf - 1 ) / e );
    vr21 = sqrt( gm_Sole / p ) * e * sin(theta21);
    vr22 = sqrt( gm_Sole / p ) * e * sin(theta22);
    if vr21 > 0 && vr22 < 0
        theta2 = theta21;
    elseif  vr21 < 0 && vr22 > 0
        theta2 = theta22;
    else 
        error('errore nella risoluzione del calcolo dell anomalia vera al primo encounter')
    end
else 
    error('non avviene trasferta poichè orbita di trasferta è circolare')
end

% velocità radiale all'encounter [km/s]
vr2 = sqrt( gm_Sole / p ) * e * sin(theta2);
% velocità tangenziale all' encounter [km/s]
vt2 = sqrt( gm_Sole / p ) * ( 1 + e * cos(theta2) );

% determinazione rx_f , ry_f , rz_f , vx_f , vy_f , vz_f all'encounter
theta_t_2 = theta2 + w; % [rad]
if i == 0 % casistica traiettoria planare
    rx_f = rf * cos(theta_t_2);
    ry_f = rf * sin(theta_t_2);
    rz_f = 0;
    vx_f = vr2 * cos(theta_t_2) - vt2 * sin(theta_t_2);
    vy_f = vr2 * sin(theta_t_2) + vt2 * cos(theta_t_2);
    vz_f = 0;
else  % casistica traiettoria 3D
    rx_f = rf * ( cos(theta_t_2)*COM1 - sin(theta_t_2)*cos(i)*SOM1 );
    ry_f = rf * ( cos(theta_t_2)*SOM1 + sin(theta_t_2)*cos(i)*COM1 );
    rz_f = rf * ( sin(theta_t_2)*sin(i) );
    vx_f = sqrt( gm_Sole / p ) * ( -COM1*( sin(theta_t_2) + e*sin(w) ) - SOM*cos(i)*( cos(theta_t_2) + e*cos(w) ) );
    vy_f = sqrt( gm_Sole / p ) * ( COM1*cos(i)*( cos(theta_t_2) + e*cos(w) ) - SOM*( sin(theta_t_2) + e*sin(w) ) );
    vz_f = sqrt( gm_Sole / p ) * sin(i)*( cos(theta_t_2) + e*cos(w) );
end




