function isit = isInteger(v)
%% Check if v is integer with tolerance tol
tol = 1e-3;
err = abs(v - round(v));
isit = (err < tol);