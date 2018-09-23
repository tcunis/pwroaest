function sp = uplus(sp)
% Unary plus operator.
%
% See UPLUS.

sp.f = cellfun(@(c) +c, sp.f, 'UniformOutput', false);

end