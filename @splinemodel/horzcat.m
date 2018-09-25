function sp = horzcat(sp1,sp2,varargin)
% Concatenates spline models vertically.

switch nargin
case 0, sp = [];
case 1, sp = sp1;

case 2
    if ~isa(sp1,'splinemodel')
        sp = sp2;
        sp.f = cellfun(@(c) horzcat(sp1,c), sp2.f, 'UniformOutput', false);
    elseif ~isa(sp2,'splinemodel')
        sp = sp1;
        sp.f = cellfun(@(c) horzcat(c,sp2), sp2.f, 'UniformOutput', false);
    else
        sp = combine(sp1,sp2,@horzcat);
    end

otherwise
    sp = horzcat(horzcat(sp1,sp2),horzcat(varargin{:}));
end