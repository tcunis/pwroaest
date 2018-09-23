function sp = subs(sp, old, new)
% See POLYNOMIAL.SUBS.

sp.H = subs(sp.H, old, new);
sp.f = cellfun(@(c) subs(c, old, new), sp.f, 'UniformOutput', false);

end