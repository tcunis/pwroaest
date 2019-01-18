function [sosc,N] = sproav_continuity(V, H, Adj, z, cstr)
% Returns SOS constraint for continuity.

if ~exist('cstr','var')
    cstr = 'c';
end

k = length(V);

I = 1:k;

soscH = cell(k,k);
soscR = cell(k,k);

I2 = false(k);
for i=I
    for j=I(Adj(i,:))
        ra = sosmdecvar(sprintf('%s%s%d%d',cstr,'a',i,j), z, length(H{i}));
        rb = sosmdecvar(sprintf('%s%s%d%d',cstr,'b',i,j), z, length(H{j}));

        soscH{i,j} = -((V{i}-V{j}) - H{i}'*ra - H{j}'*rb) >= 0;
        soscR{i,j} = [ra; rb] >= 0;
        I2(i,j) = true;
    end
end

sosc = [
    vertcat(soscH{I2})
    vertcat(soscR{I2})
];

% number of constraints
N = sum(I2(:));