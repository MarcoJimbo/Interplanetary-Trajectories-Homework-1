function dstate = Dyn_Eqs(t, state)
global gamma_max W mu_ur

% Estrazione stati
r = state(1);
% theta = state(2); % Non serve esplicitamente per le derivate
vr = state(3);
vtheta = state(4);
% C = state(5); % Non serve per le derivate

% Estrazione costati
Pr = state(6);
Ptheta = state(7);
Pvr = state(8);
Pvtheta = state(9);
PC = state(10);

Pv = [Pvr; Pvtheta];
norm_Pv = norm(Pv);

% Controllo ottimo
if (norm_Pv + PC) > 0
    % Motore Acceso
    u_r = Pvr / norm_Pv;       % Coseno direttore radiale
    u_theta = Pvtheta / norm_Pv; % Coseno direttore tangenziale
    gamma_r = gamma_max * u_r;
    gamma_theta = gamma_max * u_theta;
    dC = gamma_max;
else
    % Motore Spento
    gamma_r = 0;
    gamma_theta = 0;
    dC = 0;
end

% Equazioni del moto
dr = vr;
dtheta = vtheta/r;
dvr = vtheta^2/r - mu_ur/r^2 + gamma_r;
dvtheta = -vr*vtheta/r + gamma_theta;

% Equazioni dei costati
dPr = Ptheta*vtheta/r^2 + Pvr*(vtheta^2/r^2 + 2*mu_ur/r^3) - Pvtheta*(vr*vtheta/r^2);
dPtheta = 0;
dPvr = -Pr + Pvtheta*vtheta/r;
dPvtheta = -Ptheta/r - 2*Pvr*vtheta/r + Pvtheta*vr/r;
dPC = -(gamma_max/W) * (norm_Pv + PC);

dstate = [dr; dtheta; dvr; dvtheta; dC; dPr; dPtheta; dPvr; dPvtheta; dPC];
end

