function f = fun(X)
[~ , R] = Vincoli(X);
w = [50 , 1 , 1 , 1 , 1]';
f = sum(w.* (R.^2));
end
