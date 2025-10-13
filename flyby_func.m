function [V_plus,Delta_fin] = flyby_func(altitude,Pianeta,V_minus,Delta_in,Side_of_rotation)
%INPUT:1°altitudine:in km dalla superficie.
%2°Pianeta: Stringa ex.('Uranus') per estrarre i parametri del corpo centrale.
%(Riferirsi ai nomi dei pianeti usati nel CSPICE.)
%3°V_minus:Modulo della velocità nel sistema eliocentrico all'entrata dell'SOI.
%4°Delta_in: Angolo tra velocità del pianeta e quella dello Spacecraft all'entrata dell'SOI.
%5°Side of rotation:input 1 if Clockwise;input 2 if counterclockwise.
%OUTPUT:1°V_Plus:Modulo della velocità nel sistema eliocentrico all'uscita dell'SOI.
%2°Delta_fin:Angolo tra velocità del pianeta e quella dello Spacecraft all'uscita dell'SOI.
%!!!!!QUESTA FUNZIONE USA CSPICE, BISOGNA AVERE I KERNELS DEL PIANETA DI CUI SI VUOLE FARE IL FLYBY E AGGIORNARE mykernels.furnsh!!!!!!

% Kernels loading
cspice_furnsh( '.\kernels\mykernels.furnsh' );

%Parametri
GM_Sun = 1.32712440018e11;
GM_Planet= cspice_bodvrd(Pianeta,'GM',1);

Radius_planet=cspice_bodvrd(Pianeta, 'RADII', 3 );
Radius_planet=Radius_planet(2);%[km]
Pericenter_Radius=Radius_planet+altitude;

Data_Riferimento=cspice_str2et('2008 OCT 30 17:14:37.68');
Sun_Planet_distance=norm(cspice_spkpos(Pianeta,Data_Riferimento,'J2000','NONE','SUN'));

Planet_abs_velocicty  = sqrt(GM_Sun/Sun_Planet_distance); 
R_SoI = Sun_Planet_distance*(GM_Planet/GM_Sun)^(2/5); 

%calcolo v_inf:
v_inf=sqrt(V_minus^2+Planet_abs_velocicty^2-2*V_minus*Planet_abs_velocicty*cosd(Delta_in));
%calcolo angolo alpha-
s_a_m=(V_minus*sind(Delta_in))/v_inf;
c_a_m=(V_minus*cosd(Delta_in)-Planet_abs_velocicty)/v_inf;

alpha_minus=atan2d(s_a_m,c_a_m);

%Calcolo delta angolo tra le due v inf:
delta=asind(GM_Planet/(GM_Planet+Pericenter_Radius*v_inf^2));

%calcolo alpha+
switch Side_of_rotation
case 1 %Clock
    alpha_plus=alpha_minus+delta;
case 2 %Counterclock
    alpha_plus=alpha_minus-delta;
end

%Calcolo V_plus
V_plus=sqrt(v_inf^2+Planet_abs_velocicty^2+2*v_inf*Planet_abs_velocicty*cosd(alpha_plus));

%Calcolo Delta_fin
Delta_fin=atan2d((v_inf*sind(alpha_plus)/V_plus),((Planet_abs_velocicty+v_inf*cosd(alpha_plus))/V_plus));
end


