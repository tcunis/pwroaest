function sp = mrdivide(sp,b)
% Binary right-division operator.
%
% See MRDIVIDE.

if ~isa(b,'splinemodel')
    sp.f = cellfun(@(c) mrdivide(c,b), sp.f, 'UniformOutput', false);
else
    error('Division by spline model not supported.');
end