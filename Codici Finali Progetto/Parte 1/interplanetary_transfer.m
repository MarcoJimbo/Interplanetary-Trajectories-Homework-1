function [rx_f,ry_f,rz_f,vx_f,vy_f,vz_f,theta_t_2,delta_t,ok,status_msg,info] = interplanetary_transfer(rx_i,ry_i,rz_i,vx_i,vy_i,vz_i,gm_Sole,rf)
% Stato iniziale eliocentrico -> stato all’intersezione r=rf (PRIMO encounter)
% dedotto dallo stato iniziale (outward: vr>0; inward: vr<0).
% Funziona per ellittiche, paraboliche, iperboliche. e=0 gestito esplicitamente.
% Non solleva errori: usa ok/status_msg. Restituisce anche il tempo di volo delta_t [s].

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
%%% 8) delta_t = [s] anomalia vera + argomento al pericentro di encounter
%%% 9) ok = [bool] esito
%%% 10)status_msg = [str] diagnostica
%%% 11)info = [struct] parametri della trasferta (conic,e,p,a,rp,ra,i,theta1,theta2,branch,rf_chk,delta_t)



ok = true; status_msg = "ok";
rx_f=NaN; ry_f=NaN; rz_f=NaN; vx_f=NaN; vy_f=NaN; vz_f=NaN; theta_t_2=NaN; delta_t=NaN;

% input & grandezze di base
r_vec = [rx_i; ry_i; rz_i];
v_vec = [vx_i; vy_i; vz_i];
mu = gm_Sole;

r  = norm(r_vec);      v  = norm(v_vec);
if r==0 || ~isfinite(r) || ~isfinite(v)
    ok=false; status_msg="Input non valido (r o v non finiti)."; info=struct(); return;
end

h_vec = cross(r_vec,v_vec);  h = norm(h_vec);
if h==0
    ok=false; status_msg="Momento angolare nullo (traiettoria degenere)."; info=struct(); return;
end

r_hat = r_vec/r;
h_hat = h_vec/h;
t_hat = cross(h_hat,r_hat); t_hat = t_hat/norm(t_hat);

p   = h^2/mu;                % semilato retto
eps = v^2/2 - mu/r;          % energia specifica
a   = -mu/(2*eps);           % a=Inf per parabolica (eps≈0)

% vettore eccentricità (valido per tutte le coniche)
e_vec = (cross(v_vec,h_vec))/mu - r_hat;
e     = norm(e_vec);

% tolleranze
tol_e   = 1e-10;
tol_r   = 1e-6 * max(1,rf);
tol_den = 1e-12;

% classificazione
if e < tol_e
    conic = "circular";
elseif e < 1 - tol_e
    conic = "elliptic";
elseif abs(e-1) <= tol_e
    conic = "parabolic";
else
    conic = "hyperbolic";
end

rp = p/(1+e);  ra = NaN;
if conic == "elliptic", ra = p/(1-e); end

% inclinazione
i = acos(max(-1,min(1,h_hat(3))));

% caso e = 0 (circolare) 
if conic == "circular"
    if abs(rf - r) > tol_r
        ok=false;
        status_msg="Orbita circolare: r=costante, rf non raggiungibile.";
    else
        ok=false;
        status_msg="Orbita circolare: infinite soluzioni a r=rf; serve il tempo-di-volo/phasing.";
    end
    info = struct('conic',conic,'e',e,'p',p,'a',Inf,'rp',r,'ra',r,'i',i);
    return;
end

% raggiungibilità geometrica 
% r(θ) = p / (1 + e cosθ)  →  C = cosθ = (p/rf - 1)/e
C = (p/rf - 1)/e;
if C < -1-1e-12 || C > 1+1e-12
    ok=false; status_msg="rf fuori dal range geometrico (|cosθ|>1).";
    info = struct('conic',conic,'e',e,'p',p,'a',a,'rp',rp,'ra',ra,'i',i,'C',C);
    return;
end
C = max(-1,min(1,C));

% Iperbola: richiedi 1+e cosθ > 0 (r>0)
if conic == "hyperbolic" && (1 + e*C) <= tol_den
    ok=false; status_msg="rf non raggiungibile sull’iperbole (1+e cosθ ≤ 0).";
    info = struct('conic',conic,'e',e,'p',p,'a',a,'rp',rp,'ra',ra,'i',i,'C',C);
    return;
end

% anomalia vera iniziale θ1 (via ê e ĥ) 
e_hat = e_vec/e;
theta1 = atan2( dot(cross(e_hat, r_hat), h_hat), dot(e_hat, r_hat) );  % ∈(-π,π]

% deduzione outward/inward dallo stato iniziale 
vr0   = dot(v_vec, r_vec)/r;
v_circ = sqrt(mu/r);
if abs(vr0) > 1e-12
    want_outward = (vr0 > 0);
else
    want_outward = (v > v_circ);
end
branch = string(ifelse(want_outward,'outward','inward'));

% radici candidate per r=rf 
theta_candidates = [acos(C), -acos(C)];      % due soluzioni: ±θ
vr_candidates    = sqrt(mu/p)*e.*sin(theta_candidates);

% tieni solo la radice col segno giusto di vr (primo encounter nello stesso verso)
if want_outward
    keep = find(vr_candidates > 0);
else
    keep = find(vr_candidates < 0);
end
if isempty(keep), keep = [1 2]; end
theta_candidates = theta_candidates(keep);

% seleziona il PRIMO encounter: Δθ>0 minimo da θ1 
dtheta = arrayfun(@(th) wrapTo2Pi(th - theta1), theta_candidates);
[~,idx_min] = min(dtheta);
theta2 = theta_candidates(idx_min);

% verifica r
rf_chk = p/(1 + e*cos(theta2));
if abs(rf_chk - rf) > max(tol_r, 1e-9*rf)
    ok=false; status_msg="Intersezione r=rf non coerente (round-off eccessivo).";
    info = struct('conic',conic,'e',e,'p',p,'a',a,'rp',rp,'ra',ra,'i',i,'theta1',theta1,'theta2',theta2,'rf_chk',rf_chk);
    return;
end

% tempo di volo Δt al PRIMO encounter 
switch conic
    case "elliptic"
        % relazione θ↔E
        beta = sqrt((1 - e)/(1 + e));
        E1 = 2*atan( beta * tan(theta1/2) );
        E2 = 2*atan( beta * tan(theta2/2) );
        M1 = E1 - e*sin(E1);
        M2 = E2 - e*sin(E2);
        % assicurati ΔM>0 e minimo (primo passaggio in avanti nel tempo)
        dM = wrapTo2Pi(M2 - M1);
        n  = sqrt(mu/a^3);
        delta_t = dM / n;   % [s]

    case "parabolic"
        % Barker: t-τ = 0.5*sqrt(p^3/mu)*(D + D^3/3), D = tan(θ/2)
        D1 = tan(theta1/2);  B1 = D1 + D1^3/3;
        D2 = tan(theta2/2);  B2 = D2 + D2^3/3;
        delta_t = 0.5*sqrt(p^3/mu) * (B2 - B1); % [s]
        if delta_t < 0, delta_t = -delta_t; end % robustezza numerica

    otherwise %"hyperbolic"
        % relazione θ↔F (anomalia iperbolica): tanh(F/2)=sqrt((e-1)/(e+1))*tan(θ/2)
        gamma = sqrt((e - 1)/(e + 1));
        F1 = 2*atanh( gamma * tan(theta1/2) );
        F2 = 2*atanh( gamma * tan(theta2/2) );
        M1 = e*sinh(F1) - F1;
        M2 = e*sinh(F2) - F2;
        nH = sqrt(mu/(-a)^3);      % “mean motion” iperbolico
        delta_t = (M2 - M1) / nH;  % [s]
        if delta_t < 0, delta_t = -delta_t; end % robustezza
end

% base perifocale e trasformazione 
P_hat = e_hat;                          % verso pericentro
Q_hat = cross(h_hat, P_hat); Q_hat = Q_hat/norm(Q_hat);
R_PQW_ECI = [P_hat, Q_hat, h_hat];      % colonne

r_PQW = [rf*cos(theta2); rf*sin(theta2); 0];
v_PQW = sqrt(mu/p) * [-sin(theta2); e + cos(theta2); 0];

r_f = R_PQW_ECI * r_PQW;
v_f = R_PQW_ECI * v_PQW;

rx_f = r_f(1); ry_f = r_f(2); rz_f = r_f(3);
vx_f = v_f(1); vy_f = v_f(2); vz_f = v_f(3);

% u = θ + ω (theta_t_2) 
n_vec = cross([0;0;1], h_hat); n = norm(n_vec);
if n > 1e-12
    n_hat = n_vec/n;
    OM  = atan2(n_hat(2), n_hat(1));  %#ok<NASGU>
    w   = atan2(dot(cross(n_hat,P_hat),h_hat), dot(n_hat,P_hat));
    theta_t_2 = wrapTo2Pi(theta2 + w);
else
    theta_t_2 = wrapTo2Pi(atan2(r_f(2), r_f(1))); % quasi-equatoriale
end

info = struct( ...
    'conic',conic,'e',e,'p',p,'a',a,'rp',rp,'ra',ra,'i',i, ...
    'theta1',theta1,'theta2',theta2,'branch',branch,'rf_chk',rf_chk, ...
    'delta_t',delta_t);

end

% utilities 
function ang = wrapTo2Pi(ang)
    ang = mod(ang, 2*pi);
    if ang < 0, ang = ang + 2*pi; end
end
function out = ifelse(cond,a,b)
    if cond, out=a; else, out=b; end
end
