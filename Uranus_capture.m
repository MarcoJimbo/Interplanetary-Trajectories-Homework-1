function [d_V, delta_t, ok, status_msg, info] = Uranus_capture(v_inf, r_p, gm_Urano, r_SOI)
%%% Calcola l'impulso al pericentro per circolarizzare a r_p e il tempo
%%% di volo dalla frontiera della SOI di Urano al pericentro lungo
%%% l'iperbola di approccio.

%%% INPUT:
%%% 1)  v_inf = [km/s] velocità di eccesso iperbolico (inbound) all'infinito
%%% 2)  r_p = [km] raggio al pericentro dell'iperbola
%%% 3)  gm_Urano = [km^3/s^2] parametro gravitazionale di Urano
%%% 4)  r_SOI = [km] raggio della sfera di influenza (centro Urano)

%%% OUTPUT:
%%% 1) d_V = [km/s] impulso al pericentro per ottenere orbita circolare a r_p (valore negativo: frenata)
%%% 2) delta_t = [s] tempo dalla SOI al pericentro (ramo inbound)
%%% 3) ok = [bool] esito
%%% 4) status_msg = [str] diagnostica
%%% 5) info = [struct] parametri dell’iperbola (a,e,F_SOI, nH, v_p, ecc.)

ok = true; status_msg = "ok";
d_V = NaN; delta_t = NaN; info = struct();

mu = gm_Urano;

% controlli input minimi
if any(~isfinite([v_inf,r_p,r_SOI])) || v_inf <= 0 || r_p <= 0 || r_SOI <= 0
    ok=false; status_msg="Input non validi (v_inf>0, r_p>0, r_SOI>0)."; return;
end
if r_SOI <= r_p
    ok=false; status_msg="r_SOI deve essere maggiore di r_p."; return;
end


% Parametri iperbolici
% semiasse (magnitudine positiva) [km]
a = mu / (v_inf^2);
% eccentricità (da vincoli al pericentro)
e = 1 + (r_p * v_inf^2) / mu;   % e > 1

% velocità al pericentro [km/s]
v_p = sqrt( v_inf^2 + 2*mu/r_p );

% impulso per circolarizzare a r_p (stessa quota) [km/s]
v_circ = sqrt(mu / r_p);
d_V = v_circ - v_p;     % negativo: frenata

% anomalia iperbolica alla SOI (ramo positivo; inbound/outbound hanno |F| uguale) [rad]
arg_cosh = (r_SOI / a + 1) / e;
if arg_cosh < 1 - 1e-12
    ok=false; status_msg="Geometria incoerente: (r_SOI/a + 1)/e < 1."; 
    info = struct('a',a,'e',e,'arg_cosh',arg_cosh);
    return;
end
F_SOI = acosh( max(1, arg_cosh) );  % >= 0

% “mean motion” iperbolico (con a come magnitudine)
nH = sqrt( mu / a^3 );

% tempo di volo dalla SOI al pericentro (magnitudine)
delta_t = ( e*sinh(F_SOI) - F_SOI ) / nH;   % [s], sempre >= 0

% pacchetto info
info = struct( ...
    'a',a, 'e',e, 'F_SOI',F_SOI, 'nH',nH, ...
    'v_p',v_p, 'v_circ',v_circ, 'r_SOI',r_SOI );

end


