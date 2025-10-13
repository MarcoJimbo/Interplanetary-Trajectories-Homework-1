h_vec = linspace(200,2000,10);
dV_vec = linspace(2.90,5,10);
v_inf_mat = zeros(length(h_vec),length(dV_vec));

% NB: Use deltaV_to_vinf_new function

for i=1:length(h_vec)
    for j=1:length(dV_vec)
        v_inf_mat(i,j) = deltaV_to_vinf_new(dV_vec(j),h_vec(i));
    end
end

disp(array2table(v_inf_mat, ...
    'VariableNames',compose('dV_%.2f',dV_vec), ...
    'RowNames',compose('h_%d',h_vec)));