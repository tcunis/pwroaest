function sp = vertcat(sp1,sp2,varargin)
% Concatenates spline models vertically.

switch nargin
case 0, sp = [];
case 1, sp = sp1;

case 2
    if ~isa(sp1,'splinemodel')
        sp = vertacat(sp2,sp1);
    elseif ~isa(sp2,'splinemodel')
        sp = vertbcat(sp1,sp2);
    else
        sp = combine(sp1,sp2,@vertcat);
    end

otherwise
    sp = vertcat(vertcat(sp1,sp2),vertcat(varargin{:}));
end

end

function sp = vertacat(sp,b)
% Vertical concatenation from above, [b; sp].

    sp.f = cellfun(@(c) vertcat(b,c), sp.f, 'UniformOutput', false);
end

function sp = vertbcat(sp,b)
% Vertical concatenation from below, [sp; b].

    sp.f = cellfun(@(c) vertcat(c,b), sp.f, 'UniformOutput', false);
end