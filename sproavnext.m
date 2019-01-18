function [V0,c] = sproavnext(V, H, x, Adj, z, zi, opts)
% Computes viable Lyapunov-function for new domain.

if ~exist('opts', 'var')
    opts = sosoptions;
end
if iscell(zi)
    zi = zi{end};
end

assert(isempty(V{1}));

[V{1},c] = polydecvar('c',z);

Adj(2:end,2:end) = 0;

sosc = sproav_continuity(V,H,Adj,zi);

sosc(end+1) = V{1} >= 0;

[info,dopt] = sosopt(sosc,x,opts);

if info.feas
    V0 = subs(V{1},dopt);
    c  = subs(c,dopt);
else
    V0 = [];
    c  = [];
end
