function[i,d_V1,d_V2,t_durata] = Hohmann_Transfer(r_i,r_f,GM)
%%% function che simula Hohmann Transfer sia per trasferte verso pianeti
%%% esterni che interni. 

%%% INPUT:
%%% 1) r_i = [km] raggio circonferenza orbita del pianeta di partenza
%%% 2) r_f = [km] raggio circonferenza orbita del pianeta target
%%% 3) GM = [km^3/s^2] parametro gravitazionale corpo centrale

%%% OUTPUT:
%%% 1) i = 1 se trasferta verso pianeta esterno, 0 se trasferta verso pianeta interno
%%% 2) d_V1 = [km/s] variazione di velocità necessaria al pianeta di partenza
%%% 3) d_V2 = [km/s] variazione di velocità necessaria al pianeta target
%%% 4) t_durata = [s] durata trasferta

% controllo variabili di input
if nargin ~= 3
    error('Il numero degli input di Hohmann_Transfer deve essere 3 (r_i,r_f,GM)')
end
% differenziazione trasferta verso pianeta interno o esterno
if r_i < r_f % trasferta verso pianeta esterno
    % variazione di velocità necessaria al pianeta di partenza
    d_V1 = sqrt( 2*GM / (r_i+r_f) * r_f/r_i ) - sqrt(  GM/r_i );
    % variazione di velocità al pianeta target
    d_V2 = sqrt( GM/r_f ) - sqrt( 2*GM / (r_i+r_f) * r_i/r_f );
    i = 1;

elseif r_i > r_f % trasferta verso pianeta interno
    % variazione di velocità necessaria al pianeta di partenza
    d_V1 = sqrt( 2*GM / (r_i+r_f) * r_i/r_f ) - sqrt(  GM/r_i );
    % variazione di velocità al pianeta target
    d_V2 = sqrt( GM/r_f ) - sqrt( 2*GM / (r_i+r_f) * r_f/r_i ); 
    i = 0;
else 
    error('Pianeta di partenza e pianeta target devono avere orbite circolari di raggio diverso tra loro')
end
% durata transferta
t_durata =  pi * sqrt( ( (r_i+r_f) / 2 )^3 / GM );