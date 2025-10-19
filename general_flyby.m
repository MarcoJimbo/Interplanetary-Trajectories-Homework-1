function [v_inf2, def_angle, t_in, t_out, delta_t, ok, status_msg, info] = general_flyby(v_inf, r_p, mu_p, dV, r_SOI)

%%% function che determina la velocità di eccesso iperbolica in uscita e
%%% l'angolo di deflessione del powered flyby (manovra impulsiva tangenziale al pericentro).

%%% INPUT:
%%% 1) v_inf:        SCALAR [1x1] = velocità di eccesso iperbolico in entrata [km/s]
%%% 2) r_p:          SCALAR [1x1] = raggio al pericentro [km]
%%% 3) gm_planet:    SCALAR [1x1] = parametro gravitazionale pianeta [km^2/s^2]
%%% 4) dV:           SCALAR [1x1] = impulso erogato al pericentro [km/s]

%%% OUTPUT:
%%% 1) v_inf2:       SCALAR [1x1] = velocità di eccesso iperbolico in  uscita [km/s]
%%% 2) def_angle:    SCALAR [1x1] = angolo di rotazione velocità di eccesso iperbolico [km/s]


% ---------- init ----------
ok = true; status_msg = "ok";
v_inf2 = NaN; def_angle = NaN; t_in = NaN; t_out = NaN; delta_t = NaN;

% ---------- check ----------
if any(~isfinite([v_inf, r_p, mu_p, dV, r_SOI])) || v_inf < 0 || r_p <= 0 || mu_p <= 0 || r_SOI <= 0
    ok = false; status_msg = "Input non validi."; info = struct(); return;
end
if r_SOI <= r_p
    ok = false; status_msg = "r_SOI deve essere maggiore di r_p."; info = struct(); return;
end

% ---------- PRE-burn hyperbola ----------
eps_pre   = 0.5 * v_inf^2;                 % energia specifica (iperbolica)
vp_pre    = sqrt( v_inf^2 + 2*mu_p/r_p );
h_pre     = r_p * vp_pre;
a_pre     = -mu_p / (2*eps_pre);           % <0
aabs_pre  = -a_pre;
e_pre     = 1 + (r_p * v_inf^2) / mu_p;    % >1
nuinf_pre = acos( clamp(-1/e_pre, -1, 1) );% anomalia vera asintotica
nH_pre    = sqrt( mu_p / aabs_pre^3 );

% Tempo SOI -> pericentro sul ramo inbound (pre-burn)
F_SOI_pre = acosh( max(1, (r_SOI/aabs_pre + 1)/e_pre) );
t_in      = ( e_pre * sinh(F_SOI_pre) - F_SOI_pre ) / nH_pre;

% ---------- POST-burn conica ----------
vp_post   = vp_pre + dV;
eps_post  = 0.5*vp_post^2 - mu_p/r_p;

if eps_post <= 0
    % non iperbolica in uscita
    ok = false;
    status_msg = "Uscita non iperbolica (ε_post ≤ 0): nessun v∞_out, tempi post indefiniti.";
    v_inf2 = 0;
    e_post = NaN; nuinf_post = NaN; aabs_post = NaN; nH_post = NaN; F_SOI_post = NaN;
    % def_angle non definito tra asintoti
else
    v_inf2    = sqrt( 2*eps_post );
    a_post    = -mu_p / (2*eps_post);      % <0
    aabs_post = -a_post;
    e_post    = 1 + (r_p * v_inf2^2) / mu_p;
    nuinf_post= acos( clamp(-1/e_post, -1, 1) );
    nH_post   = sqrt( mu_p / aabs_post^3 );
    % Tempo pericentro -> SOI sul ramo outbound (post-burn)
    F_SOI_post = acosh( max(1, (r_SOI/aabs_post + 1)/e_post) );
    t_out      = ( e_post * sinh(F_SOI_post) - F_SOI_post ) / nH_post;
    delta_t      = t_in + t_out;

    % angolo tra gli asintoti (direzioni in PQW)
    u_in  = [  sin(nuinf_pre),   e_pre  + cos(nuinf_pre) ];
    u_out = [ -sin(nuinf_post),  e_post + cos(nuinf_post) ];
    cang  = dot(u_in,u_out) / (norm(u_in)*norm(u_out));
    def_angle = acos( clamp(cang, -1, 1) );
end

% ---------- info ----------
info = struct( ...
    'vp_pre', vp_pre, 'vp_post', vp_post, ...
    'eps_pre', eps_pre, 'eps_post', eps_post, ...
    'h_pre', h_pre, ...
    'aabs_pre', aabs_pre, 'e_pre', e_pre, 'nuinf_pre', nuinf_pre, 'nH_pre', nH_pre, 'F_SOI_pre', F_SOI_pre, ...
    'aabs_post', ifnan(aabs_post,NaN), 'e_post', ifnan(e_post,NaN), 'nuinf_post', ifnan(nuinf_post,NaN), ...
    'nH_post', ifnan(nH_post,NaN), 'F_SOI_post', ifnan(F_SOI_post,NaN), ...
    'r_SOI', r_SOI );

end

% -------- helper locali --------
function x = clamp(x,a,b), x = max(a, min(b, x)); end
function y = ifnan(x,alt), if isnan(x), y=alt; else, y=x; end;end