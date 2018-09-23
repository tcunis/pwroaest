function sp = uminus(sp)
% Unary minus operator.
%
% See UMINUS.

sp.f = cellfun(@(c) -c, sp.f, 'UniformOutput', false);

end