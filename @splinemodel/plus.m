function sp = plus(sp,b)
% Binary plus operator.
%
% See PLUS.

if ~isa(sp,'splinemodel')
    sp = plus(b,sp);
elseif ~isa(b,'splinemodel')
    sp.f = cellfun(@(c) plus(c,b), sp.f, 'UniformOutput', false);
else
    sp = combine(sp,b,@plus);
end