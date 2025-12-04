function dstate = Dyn_Eqs(t, state)

mu_ur = 5793939; %[km^3/s^2] parametro gravitazionale saturno
global gamma_max W
% W = 0.0032;

% Estrazione delle variabili di stato
r = state(1);
theta = state(2);
vr = state(3);
vtheta = state(4);
C = state(5);

% Estrazione delle variabili di costato
Pr = state(6);
Ptheta = state(7);
Pvr = state(8);
Pvtheta = state(9);
PC = state(10);

Pv = [Pvr; Pvtheta];

[gamma_r, gamma_theta] = gamma_control(Pv, PC);

% Equazioni di stato
dr = vr;
dtheta = vtheta/r;
dvr = vtheta^2/r - mu_ur/r^2 + gamma_r;
dvtheta = -vr*vtheta/r + gamma_theta;

if gamma_r == 0 && gamma_theta == 0
    dC = 0;
else
    dC = gamma_max;
end

% Equazioni di costato
% dPr = Pvr * (-2*mu/r^3);
% dPtheta = 0;
% dPvr = -Pr;
% dPvtheta = 0;
% dPC = -(gamma_max/W) * (norm(Pv) + PC);

dPr = (Ptheta*vtheta/r^2) + Pvr*((vtheta^2/r^2) + (2*mu_ur/r^3)) -Pvtheta*(vr*vtheta/r^2);
dPtheta = 0;
dPvr = -Pr + Pvtheta*(vtheta/r); 
dPvtheta = -Ptheta/r -Pvr*(2*vtheta/r) + Pvtheta*(vr/r); 
dPC = -(gamma_max/W) * (norm(Pv) + PC);

dstate = [dr; dtheta; dvr; dvtheta; dC;
          dPr; dPtheta; dPvr; dPvtheta; dPC];

end