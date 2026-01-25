function out = minCf(X) 
% questa funzione restituisce la funzione J da minimizzare
% il vettore X è così definito ---> X = [Pr0, Pvr0, Pvtheta0, Pc0, t_fin, Cf];
Cfinal = X(6); %%Cf è la sesta componente del vettore X
out = Cfinal;
end

