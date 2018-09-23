function sp = cleanpoly(sp, tol, deg)
% Cleans up polynomials of spline model.
%
% See POLYNOMIAL.CLEANPOLY.

sp.H = cleanpoly(sp.H,tol,deg);
sp.f = cellfun(@(c) cleanpoly(c,tol,deg), sp.f, 'UniformOutput', false);

end