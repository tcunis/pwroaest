function sp = subs(sp, old, new)
% Symbolic substitution of piecewise polynomials.
%
% Syntax:
%
%   sp3 = subs(sp2,x,p1)
%
% Replaces variables in spline model with polynomial. Number of domains in
% returned spline model remains unchanged.
%
%   sp3 = subs(p2,x,sp1)
%
% Replaces variables in polynomial with spline model. Structure of domains 
% of returned spline model remain unchanged.
%
%   sp3 = subs(sp2,x,sp1)
%
% Replaces variables in spline model with spline model. Number and
% structure of domains in returned spline model is changed to account for
% the resulting model.
%
% Inputs:
%   sp1,p1  No-by-Npts array of spline model, polynomials, or doubles.
%   sp2,p2  Nr-by-Nc spline model or polynomial array.
%   x   No-by-1 array of polynomial variables.
%
% See POLYNOMIAL.SUBS.

if isa(new,'splinemodel')
    sp = concatenate(sp,old,new);
    return
end

% else:
sp.H = subs(sp.H, old, new);
sp.f = cellfun(@(c) subs(c, old, new), sp.f, 'UniformOutput', false);

end