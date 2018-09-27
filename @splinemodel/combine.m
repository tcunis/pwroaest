function sp = combine(sp1,sp2,op)
% Combines two spline models under binary operator.
%
% Boundaries are not necessarily equal.

assert(isa(sp1,'splinemodel'), 'First argument must be a spline model.');
assert(isa(sp2,'splinemodel'), 'Second argument must be a spline model.');
assert(isa(op,'function_handle'), 'Operator must be a function handle.');

eqH = isequal(sp1.getH, sp2.getH);

if all(eqH(:))
    % spline models of equal boundaries
    sp = sp1;
    sp.f = cellfun(@(c1,c2) op(c1,c2), sp1.f, sp2.f, 'UniformOutput', false);
    
else
    % boundaries differ
    k1 = count(sp1);
    k2 = count(sp2);
    
    sp = splinemodel(k1*k2);
    
    for i1=1:k1
        for i2=1:k2
            sp.f(cb(i1,i2)) = {op(sp1.f{i1}, sp2.f{i2})};
            
            [J1,H1] = adjacent(sp1, i1);
            for j1=J1
                sp = set_adjacent(sp, cb(i1,i2), cb(j1,i2), H1(J1==j1));
            end
            
            [J2,H2] = adjacent(sp2, i2);
            for j2=J2
                sp = set_adjacent(sp, cb(i1,i2), cb(i1,j2), H2(J2==j2));
            end
        end
    end
end


function idx = cb(i1,i2)
% Returns combined index of i1 and i2.

    idx = (i1-1)*k1 + i2;
end

end