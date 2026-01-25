%% Calcolo dell'orbita effettiva su Urano dei parametri orbitali e delle manovre per arrivare all'orbita finale.
 %Vettore iniziale all'entarata dell'SOI
r_SC_Uranus_SOI_ENT_J2000=best_soluzione.Jupiter_Uranus_tr.r'; %ATTENZIONE VETTORI COLONNA
v_SC_Uranus_SOI_ENT_J2000=best_soluzione.Jupiter_Uranus_tr.v';
%Bisogna trovare Et tali che urano sia a distanza r_SOI dallo spacecraft
% d_sole_urano=0;
% et=cspice_str2et('2052 JUN 06 16:36:22.14');
% while abs(abs(r_SC_Uranus_SOI_ENT_J2000-d_sole_urano)-SOI_Urano)>100
%     et = et+1000;
%     d_sole_urano=cspice_spkpos('Uranus',et,'J2000','NONE','SUN');
% end
% disp(cspice_et2utc(et,'C',2))
% DATA utile 2147 OCT 01 13:56:22.14 per precisione di 100km oppure 2052
% JUL 06 16:36:22.14 per precisione di 1000 km
ET_SC_Uranus_SOI_ENT=cspice_str2et('2147 OCT 01 13:56:22.14');
R_J2000_BF_UR_ETentrance = cspice_pxform( 'J2000', 'IAU_URANUS', ET_SC_Uranus_SOI_ENT);

X_SC_Uranus_SOI_ENT_BFURANUS=[R_J2000_BF_UR_ETentrance*r_SC_Uranus_SOI_ENT_J2000;R_J2000_BF_UR_ETentrance*v_SC_Uranus_SOI_ENT_J2000];
 %Parametri orbitali
oe_1=cspice_oscelt(X_SC_Uranus_SOI_ENT_BFURANUS,ET_SC_Uranus_SOI_ENT,gm_Urano);

fprintf('VETTORE DI STATO E PARAMETRI ORBITALI NEL PUNTO DELLA PRIMA MANOVRA \n')
fprintf('Spacecraft Position: %e Km %e Km %e Km\n', X_SC_Uranus_SOI_ENT_BFURANUS(1:3))
fprintf('\n')
fprintf('Spacecraft Velocity: %e Km/s %e Km/s %e Km/s\n', X_SC_Uranus_SOI_ENT_BFURANUS(4:6))
fprintf('\n')
fprintf('\n')
fprintf('Perifocal distance or radius of pericenter (r_p): %f Km\n', oe_1(1))
fprintf('\n')
fprintf('Semi-major axis (a): %f Km\n', oe_1(1)/(1-oe_1(2)))
fprintf('\n')
fprintf('Eccentricity (e): %f \n', oe_1(2))
fprintf('\n')
fprintf('Inclination (i): %f deg \n', rad2deg(oe_1(3)))
fprintf('\n')
fprintf('Right Ascension of the Ascending Node (%s): %f deg\n', char(937), rad2deg(oe_1(4)))
fprintf('\n')
fprintf('Argument of Pericenter (%s): %f deg \n', char(969), rad2deg(oe_1(5)))
fprintf('\n')
fprintf('Mean Anomaly (M): %f deg\n', rad2deg(oe_1(6)))
fprintf('\n')