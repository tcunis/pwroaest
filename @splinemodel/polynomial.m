function sp = polynomial(sp)
% Creates a spline model of polynomials.

    sp.f = cellfun(@(c) polynomial(c), sp.f, 'UniformOutput', false);
    sp.H = polynomial(sp.H);
    
end