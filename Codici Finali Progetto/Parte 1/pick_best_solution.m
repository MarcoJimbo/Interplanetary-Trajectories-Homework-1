function [best, idx] = pick_best_solution(SOLUZIONE_T_M_G_U)
%PICK_BEST_SOLUTION Restituisce la soluzione con dV_tot minimo
%   [best, idx] = pick_best_solution(S) cerca dentro la cell array S
%   (tipicamente SOLUZIONE_T_M_G_U) la struct con campo .d_V_tot minimo,
%   ignorando celle vuote o con d_V_tot non valido. 
%   OUTPUT:
%     best : struct (la cella corrispondente alla soluzione migliore)
%     idx  : (opzionale) struct con indici i,j,l,k,m e dV_min
%
%   Esempio:
%     [best, idx] = pick_best_solution(SOLUZIONE_T_M_G_U);

    if nargin ~= 1
        error('pick_best_solution richiede 1 input: la cella delle soluzioni.');
    end
    S = SOLUZIONE_T_M_G_U;

    if ~iscell(S) || isempty(S)
        error('L''input deve essere una cell array non vuota.');
    end

    % Estrai dV_tot come matrice numerica (NaN per celle non valide)
    dV_grid = cellfun(@(c) get_dVtot_or_nan(c), S);

    % Trova il minimo ignorando NaN
    [min_dV, lin_idx] = min(dV_grid(:), [], 'omitnan');

    if isempty(min_dV) || isnan(min_dV)
        error('Nessuna soluzione fattibile trovata (tutte NaN).');
    end

    % Indici della cella migliore
    sz = size(S);
    [i_best, j_best, l_best, k_best, m_best] = ind2sub(sz, lin_idx);

    % Restituisci la struct corrispondente
    best = S{i_best, j_best, l_best, k_best, m_best};

    if nargout > 1
        idx = struct('i', i_best, 'j', j_best, 'l', l_best, 'k', k_best, 'm', m_best, ...
                     'dV_min', min_dV);
    end
end

% helper locale
function val = get_dVtot_or_nan(c)
    if isempty(c) || ~isstruct(c) || ~isfield(c,'d_V_tot') || ~isfinite(c.d_V_tot)
        val = NaN;
    else
        val = c.d_V_tot;
    end
end
