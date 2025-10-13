function[rx_f,ry_f,rz_f,vx_inf,vy_inf,vz_inf] = interplanetary_transfer(rx_i,ry_i,rz_i,vx_i,vy_i,vz_i,gm_Sole,rf)

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
r_versore = r_i / r;
h_versore = h_vec / h;
theta_versore = cross(h_versore,r_versore) / norm( cross(h_versore,r_versore) );
R_LVLH_ECI_t0 = [ r_versore , theta_versore , h_versore];
i = acos(h_versore(3)); % [rad] inclinazione orbita eliocentica
if i == 0
    Cot1 = R_LVLH_ECI_t0(1,1);
    Sot1 = R_LVLH_ECI_t0(2,1);
    theta_t_1 = 2 * atan( Sot1 / ( 1 + Cot1 ) );
    w = theta_t_1 - theta1; % [rad] 
    % OMEGA non definito
else
    
end

% determinazione anomalia vera al primo encounter [rad]
if vr > 0 % trasferta verso pianeta esterno
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
elseif vr< 0 % trasferta verso pianeta interno
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
end

% velocità radiale all'encounter [km/s]
vr2 = sqrt( gm_Sole / p ) * e * sin(theta2);
% velocità tangenziale all' encounter [km/s]
vt2 = sqrt( gm_Sole / p ) * ( 1 + e * cos(theta2) );
% velocità all'encounter [km/s]
vf = sqrt( vr2^2 + vt2^2 );

% determinazione matrice di rotazione da LVLH a ECI
R_LVLH_ECI = 

