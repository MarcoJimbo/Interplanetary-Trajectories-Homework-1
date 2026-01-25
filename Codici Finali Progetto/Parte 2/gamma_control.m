function [gamma_r, gamma_theta] = gamma_control(Pv, Pc)
% Questa funzione mi restituisce in output le componenti dell'accelerazione
% ottimale e nella direzione ottimale ottenute tramite la Maximum Condition
% del problema di Pontryagin
global gamma_max
% Definisco direzione ottimale della spinta
if norm(Pv) > 0
    D_opt = Pv / norm(Pv);
else
    D_opt = [0; 0];
end
if (norm(Pv) + Pc) > 0
    gamma_r = gamma_max * D_opt(1);
    gamma_theta = gamma_max * D_opt(2);
else
    gamma_r = 0;
    gamma_theta = 0;
end
end

