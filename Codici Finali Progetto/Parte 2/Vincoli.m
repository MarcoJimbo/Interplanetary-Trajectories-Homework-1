function [c, ceq] = Vincoli(X)
global gamma_max R_start r_umbriel mu_ur
% Estrazione Variabili
Pr0 = X(1); Pvr0 = X(2); Pvtheta0 = X(3); PC0 = X(4); t_fin = X(5);

% Condizioni Iniziali
r0 = R_start;

e0=0;
a0=r0/(1-e0);

theta0 = 0; 
vr0 = 0;
vtheta0 = sqrt((mu_ur*(1+e0))/r0); 

X0 = [r0, theta0, vr0, vtheta0, 0]';
P0 = [Pr0; 0; Pvr0; Pvtheta0; PC0];
IN = [X0; P0];

options = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);

try
    [~, out] = ode113(@(t, s) Dyn_Eqs(t, s), [0, t_fin], IN, options);
    
    % Dati finali
    r_f = out(end, 1);
    vr_f = out(end, 3);
    vtheta_f = out(end, 4);
    
    % Costati finali
    Pr_f = out(end, 6);
    Pvr_f = out(end, 8);
    Pvtheta_f = out(end, 9);
    PC_f = out(end, 10);
    
    % Hamiltoniana Finale
    Pv_f = [Pvr_f; Pvtheta_f];
    SF_f = norm(Pv_f) + PC_f;
    if SF_f > 0, gamma = gamma_max; else, gamma = 0; end
    
    H_f = Pr_f*vr_f + Pvr_f*(vtheta_f^2/r_f - mu_ur/r_f^2) + ...
          Pvtheta_f*(-vtheta_f*vr_f/r_f) + gamma*SF_f;
      
    % Vincoli Uguaglianza
    v_encelado_target = sqrt(mu_ur / (r_umbriel - 2000) );
    
    ceq = [
        (r_f - r_umbriel - 2000) / 1000;      % Raggio finale (scalato)
        vr_f*100;                       % Vel. radiale 0 (scalata)
        (vtheta_f - v_encelado_target)*10; % Vel. tangenziale (scalata)
        H_f;                            % H = 0
        PC_f + 1;                       % Trasversalità
    ];
    
    % Vincoli Disuguaglianza (Collisione)
    r_vector = out(:, 1);
    c = (25362 + 2000) - min(r_vector); % Se scendo sotto r_saturno, c > 0 (violation)

catch
    % Recovery robusto: se ode fallisce, restituisci penalità altissime
    ceq = ones(5,1)*1e9;
    c = 1e9; 
end
end