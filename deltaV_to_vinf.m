function [v_inf, delta_t, ok, status_msg, info] = deltaV_to_vinf(delta_V, r_p, mu_E, r_SOI_E)
% Converte un impulso tangenziale applicato in LEO in v_inf eliocentrico
% planetocentrico (rispetto alla Terra) e calcola il TEMPO di volo dal
% pericentro (punto di applicazione dell'impulso) fino alla frontiera
% della SOI terrestre lungo la traiettoria di uscita.
%
% INPUT
%   delta_V   [km/s]   impulso tangenziale su orbita LEO (>=0 per fuga)
%   h_LEO     [km]     quota LEO
%   R_E       [km]     raggio terrestre
%   mu_E      [km^3/s^2] mu_Terra
%   r_SOI_E   [km]     raggio della SOI terrestre (centrata sulla Terra)
%
% OUTPUT
%   v_inf     [km/s]   velocità di eccesso iperbolico (>=0)
%   delta_t   [s]      tempo per raggiungere la SOI (pericentro -> SOI)
%   ok        [bool]   esito
%   status_msg[string] diagnostica
%   info      [struct] parametri utili (v_c, v_p, v_esc, a, e, F_SOI, ...)
%
% NOTE
% - Se delta_V è insufficiente, ok=false e v_inf=NaN.
% - Gestisce il caso limite PARABOLICO (v_inf ~ 0) con la formula di Barker.

ok = true; status_msg = "ok";
v_inf = NaN; delta_t = NaN; info = struct();

% controlli input 
if any(~isfinite([delta_V, r_p, mu_E, r_SOI_E])) || ...
        (r_p <= 0) || (mu_E <= 0) || (r_SOI_E <= 0)
    ok=false; status_msg="Input non validi."; return;
end
if r_SOI_E <= (r_p)
    ok=false; status_msg="r_SOI_E deve essere > r_LEO."; return;
end

% grandezze di base  
v_c   = sqrt(mu_E / r_p);      % [km/s] velocità circolare in LEO
v_p   = v_c + delta_V;         % [km/s] velocità dopo l'impulso
v_esc = sqrt(2*mu_E / r_p);    % [km/s] velocità di fuga parabolica

% caso impulso negativo (non fuga): non bloccare, ma segnala
if delta_V < 0
    ok=false; status_msg="delta_V negativo: traiettoria non di fuga.";
    info = struct('r_p',r_p,'v_c',v_c,'v_p',v_p,'v_esc',v_esc);
    return;
end

% classificazione: parabolica vs iperbolica 
tol = 1e-8; % [km/s] tolleranza su v_inf
if v_p < v_esc - tol
    ok=false; status_msg="Delta-V insufficiente: v_p < v_esc (no fuga).";
    info = struct('r_p',r_p,'v_c',v_c,'v_p',v_p,'v_esc',v_esc);
    return;
end

% v_inf (se iperbolica) o ~0 (parabolica)
v_inf_raw = sqrt(max(0, v_p^2 - v_esc^2));  % = sqrt(v_p^2 - 2 mu/r_p)
is_parabolic = (v_p <= v_esc + tol) && (v_inf_raw < 10*tol);
v_inf = v_inf_raw;

% tempo per raggiungere la SOI 
if is_parabolic
    % PARABOLA di fuga (Barker) 
    % h = r_p * v_esc (pericentro)
    h = r_p * v_esc;
    p = h^2 / mu_E;
    % r(θ) = p / (1 + cosθ)  =>  cosθ_SOI = p/r_SOI - 1  (clamp)
    C = p/r_SOI_E - 1;
    C = max(-1, min(1, C));
    theta_SOI = acos(C);
    % tempo da pericentro: t = 0.5*sqrt(p^3/mu)*(D + D^3/3), D=tan(θ/2)
    D  = tan(theta_SOI/2);
    B  = D + D^3/3;
    delta_t = 0.5 * sqrt(p^3/mu_E) * B;  % [s]
    branch = "parabolic";

    a = Inf; e = 1; F_SOI = NaN; nH = 0;

else
    % IPERBOLA di fuga 
    % parametri iperbolici
    a = mu_E / v_inf^2;                 % magnitudine semiasse (>0)
    e = 1 + (r_p * v_inf^2) / mu_E;     % eccentricità > 1

    % anomalia iperbolica alla SOI (r = a(e coshF - 1))
    arg_cosh = (r_SOI_E / a + 1) / e;
    if arg_cosh < 1 - 1e-12
        ok=false; status_msg="Geometria incoerente: (r_SOI/a + 1)/e < 1.";
        info = struct('r_p',r_p,'v_c',v_c,'v_p',v_p,'v_esc',v_esc,'a',a,'e',e,'arg_cosh',arg_cosh);
        return;
    end
    F_SOI = acosh(max(1, arg_cosh));
    nH    = sqrt(mu_E / a^3);

    % tempo dal pericentro alla SOI
    delta_t = ( e*sinh(F_SOI) - F_SOI ) / nH;   % [s]
    branch = "hyperbolic";
end

% info di ritorno 
info = struct( ...
    'branch',  branch, ...
    'r_p',     r_p, ...
    'v_c',     v_c, ...
    'v_p',     v_p, ...
    'v_esc',   v_esc, ...
    'a',       a, ...
    'e',       e, ...
    'F_SOI',   ifelse(isfinite(v_inf) && ~is_parabolic, F_SOI, NaN), ...
    'nH',      ifelse(isfinite(v_inf) && ~is_parabolic, nH, 0), ...
    'r_SOI',   r_SOI_E );

end

% utility locale 
function out = ifelse(cond, a, b)
    if cond, out = a; else, out = b; end
end
