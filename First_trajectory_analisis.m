clear all
clc

% Kernels:
cspice_furnsh( '.\kernels\mykernels.furnsh' );

%Anlalisi Homan per vincoli su velocità iniziali
Data_Riferimento=cspice_str2et('2009 APR 27 17:14:37.68');
gm_sun=cspice_bodvrd('Sun','GM',1);
gm_earth=cspice_bodvrd('Earth','GM',1);

h_LEO=300;
r_LEO= 1.495978707e8+h_LEO; %km

vesc=sqrt(2*gm_earth/r_LEO);

r_p_earth=149597870; %km

r_p_mars= 227936640 ; %norm(cspice_spkpos('Mars',Data_Riferimento,'J2000','NONE','SUN'));

r_p_jup=778412010; %norm(cspice_spkpos('Jupiter',Data_Riferimento,'J2000','NONE','SUN'));

r_p_ura=2870972200;

a_h_1=(r_p_mars+r_p_earth)/2;

a_h_2=(r_p_earth+r_p_jup)/2;

a_h_3=(r_p_earth+r_p_ura)/2;

V_per_earth_mars=sqrt(gm_sun*((2/r_p_earth)-(1/a_h_1)));
V_per_earth_jup=sqrt(gm_sun*((2/r_p_earth)-(1/a_h_2)));
V_per_earth_ura=sqrt(gm_sun*((2/r_p_earth)-(1/a_h_3)));

VP1=sqrt(gm_sun/r_p_earth);

Dv_min=V_per_earth_mars-VP1; % Velocità che nel sistema terra sarebbe vinf da aggiungere a quella della terra per avere la velocità finale nel sistema eliocentrico.

DELTA_DA_LEO_min=sqrt(Dv_min^2+vesc^2)-sqrt(gm_earth/r_LEO); %2.894294522547054 nello scrpt di dario

Dv_max=V_per_earth_jup-VP1;
vP_max=sqrt(Dv_max^2+vesc^2);
DELTA_DA_LEO_max=vP_max-sqrt(gm_earth/r_LEO);

Dv_ura=V_per_earth_ura-VP1;
%Vettore Delta v:
veradv_min=2.894294522547054; %Per lo script di Dario il Dv minimo è leggermente maggiore di quello teorico.
Vettore_deltaV=linspace(veradv_min,DELTA_DA_LEO_max,100);

Mars_arrvival_cond=zeros(100,5);
%La matrice esprime nei 100 dv analizzati le condizioni di arrivo su marte.
for i=1:length(Vettore_deltaV)
    [V_inf_terra,V_ext_terra_SoI,True_anomaly_Mars_arrival,V_minus,Delta_i_ing]=arco_1(Vettore_deltaV(i),h_LEO);
    Mars_arrvival_cond(i,:)=[V_inf_terra,V_ext_terra_SoI,True_anomaly_Mars_arrival,V_minus,Delta_i_ing];
end
%!!!!!!CAPIRE SE L'ANGOLO DELTA I di entrata ha senso o meno, perchè non
%sono sicuro che con atan2 lo calcoli bene.

%Calcolare i flyby per ogni DeltaV e dati di uscita dalla soi di marte ,
%clacolare la conica fino alla Soi di Giove e 








                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               