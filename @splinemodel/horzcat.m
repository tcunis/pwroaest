function sp = horzcat(sp1,sp2,varargin)
% Concatenates spline models vertically.

switch nargin
case 0, sp = [];
case 1, sp = sp1;

case 2
    if ~isa(sp1,'splinemodel')
        sp = horzlcat(sp2,sp1);
    elseif ~isa(sp2,'splinemodel')
        sp = horzrcat(sp1,sp2);
    else
        sp = combine(sp1,sp2,@horzcat);
    end

otherwise
    sp = horzcat(horzcat(sp1,sp2),horzcat(varargin{:}));
end

end

function sp = horzlcat(sp,b)
% Vertical concatenation from left, [b sp].

    sp.f = cellfun(@(c) horzcat(b,c), sp.f, 'UniformOutput', false);
end

function sp = horzrcat(sp,b)
% Vertical concatenation from right, [sp b].

    sp.f = cellfun(@(c) horzcat(c,b), sp.f, 'UniformOutput', false);
end