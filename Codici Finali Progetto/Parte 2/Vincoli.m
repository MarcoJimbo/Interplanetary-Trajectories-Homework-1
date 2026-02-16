function [c,ceq] = Vincoli(X)
global R_start r_umbriel mu_ur W

L  = R_start;
V  = sqrt(mu_ur *(1+0.6)/ L);
T  = L / V;

% Decision variables: [P1_0bar, P3_0bar, P4_0bar, P5_0, tfbar]
P1_0 = X(1);
P3_0 = X(2);
P4_0 = X(3);
P5_0 = X(4);
tfbar = X(5)/T;               % <-- tf in tempo normalizzato

% IC bar (se L = R_start)
r0  = 1;
th0 = 0;
vr0 = 0;
vt0 = 1;                    % perché V = sqrt(mu/L)
x50 = 1;

IN = [r0; th0; vr0; vt0; x50; P1_0; 0; P3_0; P4_0; P5_0];


opts = odeset('RelTol',1e-9,'AbsTol',1e-11,'NormControl','on','Events',@guard); % meglio senza Events in ottimizzazione

warnId = 'MATLAB:ode113:IntegrationTolNotMet';
ws = warning('query', warnId);
warning('error', warnId);

try
    [tbar,out] = ode113(@Dyn_Eqs, [0 tfbar], IN, opts);
    if tbar(end) < 0.99 * tfbar
        pen = (1 - tbar(end)/tfbar);
        ceq = 1e3*pen*ones(5,1);
        c = [];
        return
    end
catch
    ceq = 1e4*ones(5,1);
    c = [];
    return
end
rf  = out(end,1);
vrf = out(end,3);
vtf = out(end,4);
P5f = out(end,10);

% Target bar
rt  = r_umbriel / L;
vt_target = sqrt(1/rt);     % v_circ bar con mu_bar=1

% Hamiltoniana bar al finale: Hbar = Pbar' * fbar(states)
xf = out(end,:).';
fbar = Dyn_Eqs(tbar(end), xf);
Pbar = xf(6:10);
Hbar = Pbar.' * fbar(1:5);

ceq = [
    rf - rt;
    vrf;
    vtf - vt_target;
    P5f + 1;    % <-- fuel-opt Mayer fondamentale
    Hbar;       % <-- solo se tfbar libero e non saturo a bound
];

c = [];
end
