function sp = mtimes(sp1,sp2)
% Binary mtimes operator.
%
% See MTIMES.

if ~isa(sp1,'splinemodel')
    sp = mltimes(sp2,sp1);
elseif ~isa(sp2,'splinemodel')
    sp = mrtimes(sp1,sp2);
else
    sp = combine(sp1,sp2,@mtimes);
end

end

function sp = mrtimes(sp,b)
%MRTIMES Right-hand side multiplication sp*b.
    sp.f = cellfun(@(c) mtimes(c,b), sp.f, 'UniformOutput', false);
end

function sp = mltimes(sp,b)
%MLTIMES Left-hand side multiplication b*sp
    sp.f = cellfun(@(c) mtimes(b,c), sp.f, 'UniformOutput', false);
end