function[dV] = dV_setter(r1,r2,gm_Sole,GM,r_p,V_p,k1,k12,dV_step)
%%% function che determina il vettore riga dV variabile di progetto.
%%% come dV_min viene scelto:
%%% 1) se trasferta verso pianeta esterno impulso che comporta trasferta alla Hohmann 
%%% (se il livello energetico non è già sufficiente dopo flyby). Se invece
%%% il livello energetico è gia sufficiente dopo flyby dV_min viene scelto
%%% nullo e dV_max = 5 [km/s]
%%% 2) se trasferta verso pianeta interno come dv_max impulso che comporta trasferta alla Hohmann 
%%% (se il livello energetico non è già sufficiente dopo flyby). Se invece
%%% il livello energetico è gia sufficiente dopo flyby dV_min viene scelto
%%% nullo e dV_max = -5 [km/s]

%%% INPUT:
%%% 1) r1 = [km] raggio orbita circolare di partenza
%%% 2) r2 = [km] raggio orbita circolare target
%%% 3) gm_Sole = [km^3/s^2] parametro gravitazionale Sole
%%% 4) GM = [km^3/s^2] parametro gravitazionale pianeta 
%%% 5) r_p = [km] raggio al pericentro di orbita di partenza intorno pianeta
%%% 6) V_p = [km/s] velocità al pericentro di orbita di partenza intorno pianeta (precedentemente manovra)
%%% 7) k1 = costante moltiplicativa che definisce dV_max = k1 * dV_min
%%% 8) k1 = [km/s] impulso massimo erogabile se con dV = 0 si raggiunge prossimo pianeta di missione
%%% 9) dV_step = [km/s] risoluzione di variabile di progetto dV

%%% OUTPUT:
%%% dV = vector 1xn [km/s] degli impulsi erogati al pericentro di orbita dentro SOI

if nargin ~= 9
    error('Il numero di input di dV_setter deve essere 9 (r1,r2,gm_Sole,GM,r_p,V_p,k1,k12,dV_step)')
end

[i,d_V1] = Hohmann_Transfer(r1,r2,gm_Sole);
if i
   dV_min = sqrt( 2*( GM / (r_p) + d_V1^2 / 2) ) - V_p ;
    if dV_min > 0
        dV_max = k1 * dV_min;
    else 
        dV_min = 0;
        dV_max = k12;
    end
else 
    dV_max = sqrt( 2*( GM / (r_p) + d_V1^2 / 2) ) - V_p ;
    if dV_max < 0
        dV_min = k1 * dV_max;
    else 
        dV_max = 0;
        dV_min = -k12;
    end
end

dV = dV_min : dV_step : dV_max;