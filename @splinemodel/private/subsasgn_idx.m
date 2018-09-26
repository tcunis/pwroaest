function sp = subsasgn_idx(sp,s,b)
% Subscripted index assignment SP(s) = A.
%

if ~isa(b,'splinemodel')
    sp.f = cellfun(@(c) csubsasgn(c,s,b), sp.f, 'UniformOutput', false);
else
    sp = combine(sp,b,@(c1,c2) csubsasgn(c1,s,c2));
end

end

function out = csubsasgn(c,varargin)
%CSUBSASIGN Wrapper function for SUBSASGN.
%
% This avoids an error if built-in subsasgn is used in anonymous function.

    out = subsasgn(c,varargin{:})
end