function sp = mldivide(b,sp)
% Binary left-division operator.
%
% See MLDIVIDE.

if ~isa(b,'splinemodel')
    sp.f = cellfun(@(c) mldivide(b,c), sp.f, 'UniformOutput', false);
else
    error('Division by spline model not supported.');
end