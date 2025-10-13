clear all
clc 
format long

h_leo_vec=(450:10:2000); %low earth orbit height limits ogni 10 km
Deltav_vec=(2.9:0.1:4.5); %[km/s] ogni 100m/s

for i=1:length(Deltav_vec)
    for j=1:length(h_leo_vec)
        v_inf(i,j)=deltaV_to_vinf_new(Deltav_vec(i),h_leo_vec(j));
    end
end

surf(h_leo_vec,Deltav_vec,v_inf)