function sp = subsref_idx(sp, s)
% Subscripted index reference A = SP(s).
%
% See also SPLINEMODE.SUBSREF, SUBSREF.

sp.f = cellfun(@(c) subsref(c,s), sp.f, 'UniformOutput', false);

end