function sp = mtimes(sp,b)
% Binary mtimes operator.
%
% See MTIMES.

if ~isa(sp,'splinemodel')
    sp = mtimes(b,sp);
elseif ~isa(b,'splinemodel')
    sp.f = cellfun(@(c) mtimes(c,b), sp.f, 'UniformOutput', false);
else
    sp = combine(sp,b,@mtimes);
end