function dstate = Dyn_Eqs(tbar, state)
global gamma_max mu_ur R_start W

L  = R_start;
V  = sqrt(mu_ur*(1+0.6) / L);
T  = L / V;
A  = V / T;                 % = mu/L^2
gbar = gamma_max / A;       % gamma_max normalizzata
kVW  = V / W;               % fattore V/W

% Stati bar
r   = state(1);
vr  = state(3);
vt  = state(4);
x5  = state(5);

% Costati bar
P1 = state(6);
P2 = state(7);
P3 = state(8);
P4 = state(9);
P5 = state(10);

% Switching (equivalente in segno)
Sigma = hypot(P3,P4)/x5 + kVW*P5;

% Controllo bar
u1 = (Sigma > 0) * gbar;

den = hypot(P3,P4);
if den < 1e-12
    su2 = 0;  cu2 = -1;
else
    su2 = -P3/den;
    cu2 = -P4/den;
end

% Dinamica bar (mu_bar = 1)
dr  = vr;
dth = vt/r;
dvr = vt^2/r - 1/r^2 + (u1/x5)*su2;
dvt = -vr*vt/r + (u1/x5)*cu2;
dx5 = -kVW * u1;            % <-- QUI il fattore V/W è fondamentale

% Costati bar (stessa struttura, mu_bar = 1)
dP1 = P2*vt/r^2 - P3*(-vt^2/r^2 + 2/r^3) - P4*(vr*vt/r^2);
dP2 = 0;
dP3 = -P1 + P4*vt/r;
dP4 = -P2/r - 2*P3*vt/r + P4*vr/r;
dP5 = -u1/x5^2 * hypot(P3,P4);

dstate = [dr; dth; dvr; dvt; dx5; dP1; dP2; dP3; dP4; dP5];
end
