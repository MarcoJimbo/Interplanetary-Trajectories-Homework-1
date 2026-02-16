function [value,isterminal,direction] = guard(t,y)
global R_start
L = R_start;
r  = y(1);
value = r - 1e5/L;   % esempi di soglie
isterminal = 1;
direction = 0;
end
