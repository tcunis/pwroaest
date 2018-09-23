function tf = eq(sp1,sp2)
% Determines whether two spline models are identical.
%
% See EQ.


if ~isa(sp1,'splinemodel') || ~isa(sp2,'splinemodel')
    tf = false;
    return
end

% else
eqf = cell2mat(cellfun(@(c1,c2) isequal(c1,c2), sp1.f, sp2.f, 'UniformOutput',false));
eqH = isequal(sp1.getH, sp2.getH);

tf = all([eqf(:); eqH(:)]);

end