function sp = power(sp,b)
% Binary power operator.
%
% See POWER.

if ~isa(b,'splinemodel')
    sp.f = cellfun(@(c) power(c,b), sp.f, 'UniformOutput', false);
else
    error('Power of spline model not supported.');
end