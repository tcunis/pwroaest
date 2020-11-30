function sp = concatenate(sp2, old, sp1)
% Concatenate spline models.
%
% See SUBS.

if isa(sp2,'polynomial')
    sp = sp1;
    sp.f = cellfun(@(c) subs(sp2, old, c), sp.f, 'UniformOutput', false);
    return
end

% else:
k1 = count(sp1);
k2 = count(sp2);

sp = splinemodel(k1*k2);

cb = reshape(1:count(sp),k1,k2);

for i1=1:k1
    for i2=1:k2
        f1 = sp1.f{i1};
        f2 = sp2.f{i2};
        sp.f(cb(i1,i2)) = {subs(f2,old,f1)};

        [J1,H1] = adjacent(sp1, i1);
        for j1=J1
            sp = set_adjacent(sp, cb(i1,i2), cb(j1,i2), H1(J1==j1));
        end

        [J2,H2] = adjacent(sp2, i2);
        H2 = subs(H2, old, f1);
        for j2=J2
            sp = set_adjacent(sp, cb(i1,i2), cb(i1,j2), H2(J2==j2));
        end
    end
end

end