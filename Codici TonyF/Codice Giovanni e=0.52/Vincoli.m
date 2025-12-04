function [c, ceq] = Vincoli(X)  % Nomi standard per fmincon: c (inequality), ceq (equality)
format long 
r_umbriel = 266000; %[km] raggio orbita encelado attorno a saturno
r_ur = 51118; %[km] raggio di saturno
mu_ur = 5793939; %[km^3/s^2] parametro gravitazionale saturno

R_start = 2*r_ur; %[km] %generalmente tra 1.000 km e 100.000 km sopra la superficie del pianeta, noi decidiamo di arrivare all'esterno della fascia F degli anelli di Saturno che hanno un raggio di 140180km dal centro di Saturno
e0=0.52;
a0=R_start/(1-e0);
ra0=(1+e0)*a0;

% vincoli iniziali
r0 = ra0;
theta0 = pi; 
vr0 = 0;
vtheta0 = sqrt(mu_ur/a0*(1-e0)/(1+e0)); 

% vincoli finali
rf = r_umbriel ; 
ef=(ra0-rf)/(ra0+rf);
af=rf/(1-ef);
vrf = 0;
vthetaf = sqrt(mu_ur/af*(1+ef)/(1-ef));  %thetaf è unparametro libero

% parametri motore
global gamma_max

% ODE options
Tol0 = 3e-14;
Tol1 = 1e-14;
options = odeset('RelTol', Tol0, 'AbsTol', Tol1);

% definisco variabili in ingresso
Pr0 = X(1);
Ptheta0 = 0;
Pvr0 = X(2);
Pvtheta0 = X(3);
PC0 = X(4);
t_fin = X(5);
Cf = X(6);

C0 = 0;
PCf = -1; % valore ricavato da eq. di trasversalità

X0 = [r0, theta0, vr0, vtheta0, C0]';
P0 = [Pr0; Ptheta0; Pvr0; Pvtheta0; PC0];
IN = [X0; P0];

try
    [t, integrator_output] = ode113(@(t, state) Dyn_Eqs(t, state), [0, t_fin], IN, options);

    % definisco tutti gli output della dynamical al tempo finale
    r_f = integrator_output(end, 1);
    theta_f = integrator_output(end, 2);
    vr_f = integrator_output(end, 3);
    vtheta_f = integrator_output(end, 4);
    C_f = integrator_output(end, 5);
    Pr_f = integrator_output(end, 6);
    Ptheta_f = integrator_output(end, 7); 
    Pvr_f = integrator_output(end, 8); 
    Pvtheta_f = integrator_output(end, 9); 
    PC_f = integrator_output(end, 10);

    % definisco moltiplicatore velocità al tempo finale
    Pv_f = [Pvr_f; Pvtheta_f];

    % controllo ottimo
    if (norm(Pv_f) + PC_f) > 0
        sgn_fun = 1;
    else
        sgn_fun = 0;
    end

    % calcolo hamiltoniana ottima al tempo finale
    Hf = Pr_f * vr_f + Pvr_f * ((vtheta_f^2/r_f) - (mu_ur/r_f^2)) + ...
         Pvtheta_f * (-vtheta_f * vr_f / r_f) + gamma_max * sgn_fun * (norm(Pv_f) + PC_f);

    % VINCOLI DI UGUAGLIANZA (devono essere = 0)
    ceq = [r_f - rf; 
          (vr_f - vrf)*100; 
          (vtheta_f - vthetaf)*10;
          Hf;           % Hamiltoniana finale = 0
          PC_f - PCf]  % Costato finale = valore desiderato

    % VINCOLI DI DISUGUAGLIANZA (devono essere <= 0)
    c = [];  % Nessun vincolo di disuguaglianza

catch ME
    % Se l'integrazione fallisce, restituisci vincoli violati
    fprintf('Errore durante l''integrazione: %s\n', ME.message);
    ceq = ones(5, 1) * 1e6;  % Valore grande per violazione
    c = [];
end

end