function [t_min, T_repeat, details] = common_wait_solver(tw, S, epsTol, dt_user,opts)
% COMMON_WAIT_SOLVER
% Trova la minima soluzione positiva t per:
%   t = tw(i) + k_i * S(i), i=1..3
% con tolleranza epsTol (in secondi), usando:
%   1) CRT esatto su interi con prova multipla di Δt
%   2) CRT "tollerante" sulla terza congruenza
%   3) Fallback: scansione pairwise (12) vs (3) entro un orizzonte
%
% INPUT
%   tw      [3x1] double   (secondi)  - t_wait iniziali
%   S       [3x1] double >0 (secondi) - periodi sinodici
%   epsTol  scalar >0       (secondi) - tolleranza
%   dt_user (opzionale)     (secondi) - Δt preferito; se NaN, scelto automaticamente
%   opts    (opzionale struct):
%           .horizon_years (default 5e6) orizzonte per fallback scan (anni)
%           .max_trials_dt (default 12)  tentativi diversi di Δt
%           .scan_repeats  (default 2e5) ripetizioni max nella scansione pairwise
%
% OUTPUT
%   t_min     : prima soluzione positiva trovata (secondi) (NaN se non trovata)
%   T_repeat  : periodo di ripetizione stimato (secondi) per la soluzione trovata
%   details   : struct con info (dt usato, metodo, residui, ecc.)

    arguments
        tw (3,1) double
        S  (3,1) double {mustBePositive}
        epsTol (1,1) double {mustBePositive}
        dt_user (1,1) double {mustBePositive} = NaN
        opts.horizon_years (1,1) double {mustBePositive} = 5e6
        opts.max_trials_dt (1,1) double {mustBeInteger, mustBePositive} = 12
        opts.scan_repeats (1,1) double {mustBeInteger, mustBePositive} = 2e5
    end

    tw = mod(tw, S);  % normalizza

    % ---- costruisci lista Δt da provare ----
    dt_list = build_dt_list(S, epsTol, dt_user, opts.max_trials_dt);

    last_err = [];
    for idt = 1:numel(dt_list)
        dt = dt_list(idt);
        N  = max(1, round(S./dt));
        W  = round(tw./dt);
        Se = N.*dt;

        % 1) CRT esatto: ((1,2) -> 12), poi con (3)
        [x12, m12, ok12] = crt_pair(W(1), N(1), W(2), N(2));
        if ~ok12
            last_err = "CRT 12 incompatibile";
            continue
        end
        [x123, m123, ok123] = crt_pair(x12, m12, W(3), N(3));
        if ok123
            % soluzione esatta su reticolo
            x = mod(x123, m123);
            [t_min, T_repeat, residuals] = finalize_time(x, m123, dt, tw, S, Se);
            if all(residuals <= epsTol + 10*eps)
                details = pack_details('CRT exact', dt, N, W, Se, residuals, m12, m123, 0, 0, true, []);
                return
            end
            % se residui > eps, prova refinement più avanti (o prova altro dt)
            last_err = "CRT exact residui > eps";
        else
            last_err = "CRT 123 incompatibile";

            % 2) CRT tollerante: allinea 12 con 3 entro Ktol tick
            g = gcd(m12, N(3));
            r = mod(W(3) - x12, g); r = min(r, g-r);           % distanza in tick
            Ktol = ceil(epsTol / dt);
            if r <= Ktol
                delta = mod(W(3) - x12, g);
                if delta > g/2, delta = delta - g; end
                x_adj = x12 + delta;
                m123  = lcm(m12, N(3));
                x     = mod(x_adj, m123);
                [t_min, T_repeat, residuals] = finalize_time(x, m123, dt, tw, S, Se);
                if all(residuals <= epsTol + 10*eps)
                    details = pack_details('CRT tolerant', dt, N, W, Se, residuals, m12, m123, Ktol, delta, true, []);
                    return
                else
                    last_err = "CRT tolerant residui > eps";
                end
            end
        end
    end

    % 3) Fallback: pairwise scan 12 vs 3 entro orizzonte
    %    (utile quando le classi residue sono "quasi" compatibili ma non perfette)
    dt = dt_list(min(2,numel(dt_list)));   % usa un dt ragionevole tra i migliori
    N  = max(1, round(S./dt));  W = round(tw./dt);  Se = N.*dt;

    [x12, m12, ok12] = crt_pair(W(1), N(1), W(2), N(2));
    if ~ok12
        t_min = NaN; T_repeat = NaN;
        details = pack_details('pairwise-scan: no 12', dt, N, W, Se, NaN(3,1), NaN, NaN, NaN, NaN, false, last_err);
        return
    end

    T12 = m12*dt;   t0 = mod(x12, m12)*dt;  if t0==0, t0=T12; end
    horizon = opts.horizon_years * 365.25 * 86400;
    Kmax = min(opts.scan_repeats, ceil(horizon / T12));

    best_res = inf; best_t = NaN; best_k = NaN;
    for k = 0:Kmax
        tc = t0 + k*T12;
        % allinea al terzo periodo reale
        k3 = round( (tc - tw(3)) / S(3) );
        t3 = tw(3) + k3*S(3);
        res3 = min( mod(tc - t3, S(3)), mod(t3 - tc, S(3)) );
        if res3 < best_res
            best_res = res3; best_t = tc; best_k = k3;
            if best_res <= epsTol, break; end
        end
    end

    if best_res <= epsTol
        t_min = best_t;
        T_repeat = lcm(m12, N(3))*dt;  % stima periodo complessivo
        residuals = [ ...
            min(mod(t_min - tw(1), S(1)), mod(tw(1) - t_min, S(1))); ...
            min(mod(t_min - tw(2), S(2)), mod(tw(2) - t_min, S(2))); ...
            min(mod(t_min - tw(3), S(3)), mod(tw(3) - t_min, S(3)))  ];
        details = pack_details('pairwise-scan', dt, N, W, Se, residuals, m12, lcm(m12,N(3)), 0, 0, true, last_err);
        return
    else
        t_min = NaN; T_repeat = NaN;
        details = pack_details('pairwise-scan: not found', dt, N, W, Se, [NaN;NaN;NaN], m12, lcm(m12,N(3)), 0, 0, false, ...
                               struct('last_err',last_err,'best_residual_sec',best_res,'best_t_sec',best_t,'best_k3',best_k,'Kmax',Kmax));
        % non trovato nell’orizzonte: aumenta epsTol, aumenta horizon_years o riduci dt
        return
    end
end

%% ---------------- helpers ----------------
function dt_list = build_dt_list(S, epsTol, dt_user, max_trials)
    if isnan(dt_user)
        base = epsTol/5;
    else
        base = dt_user;
    end
    % costruisci una piccola scala di dt attorno a "base" (più grossolano e più fine)
    scal = [1, 0.75, 0.5, 0.33, 0.25, 0.2, 0.15, 0.1, 1.25, 1.5, 2.0, 3.0];
    scal = scal(1:min(numel(scal),max_trials));
    dt_list = base * scal;

    % filtra dt troppo fini/grossolani
    dt_min = epsTol/20;            % non più fine di eps/20 (puoi regolare)
    dt_max = min(S)/1e4;           % evita moduli giganteschi
    dt_list = dt_list(dt_list>=dt_min & dt_list<=dt_max);

    % ordina per vicinanza al base
    [~,ix] = sort(abs(dt_list - base)); dt_list = dt_list(ix);
    if isempty(dt_list), dt_list = base; end
end

function [t_min, T_repeat, residuals] = finalize_time(x, m, dt, tw, S, Se)
    % converte soluzione intera -> tempo, fa micro-refinement continuo
    if x==0, t0 = m*dt; else, t0 = x*dt; end
    k  = round((t0 - tw)./S);
    ti = tw + k.*S;
    w  = 1./Se;
    t_ref = sum(w.*ti)/sum(w);
    if t_ref <= 0
        t_ref = t_ref + ceil(-t_ref/(m*dt))*m*dt;
    end
    t_min = t_ref;
    T_repeat = m*dt;
    residuals = arrayfun(@(ti,Si) min(mod(t_min-ti,Si), mod(ti-t_min,Si)), tw, S);
end

function details = pack_details(method, dt, N, W, Se, residuals, m12, m123, Ktol, delta, success, extra)
    details = struct('method',method,'dt',dt,'N',N,'W',W,'Seff',Se, ...
                     'residuals',residuals,'m12',m12,'m123',m123, ...
                     'Ktol_tick',Ktol,'shift_tick',delta,'success',success,'extra',extra);
end

function [x, m, ok] = crt_pair(a1, m1, a2, m2)
    a1 = mod(a1, m1); a2 = mod(a2, m2);
    [g, s, ~] = egcd(m1, m2);
    if mod(a2 - a1, g) ~= 0
        x = NaN; m = NaN; ok = false; return;
    end
    l = (m1/g)*m2;
    k = ((a2 - a1)/g) * s;
    k = mod(k, m2/g);
    x = mod(a1 + m1*k, l);
    x = round(x); m = l; ok = true;
end

function [g, x, y] = egcd(a, b)
    if b == 0, g = a; x = 1; y = 0;
    else
        [g1, x1, y1] = egcd(b, mod(a,b));
        x = y1; y = x1 - floor(a/b)*y1; g = g1;
    end
end
