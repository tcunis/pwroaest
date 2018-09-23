function [gmin,x] = spdistance(Hi,V,Hj,opts)

% polynomial variables
x = unique([Hi.varname; V.varname; Hj.varname]);

% default options
if ~exist('opts','var') || isempty(opts)
    opts = optimoptions('fmincon');
end

g = pvar('g');
y = [g;x];

ceq = [
    V - g
%     prod(Hi)
];

% if length(Hi) > 1
    c = [
        Hi
        Hj
    ];
% else
%     c = Hj;
% end

problem.solver    = 'fmincon';
problem.objective = @(z) z(1);
problem.x0        = zeros(size(y));
problem.lb        = [0; -Inf; -Inf];
problem.nonlcon  = @(z) spdistcon(ceq,c,y,z);
problem.options   = opts;

[y,gmin,info] = fmincon(problem);

if info == 1
    x = y(2:end);
else
    x    = [];
    gmin = [];
end

end


function [c,ceq] = spdistcon(ceq,c,y,z)
    ceq = double(subs(ceq,y,z));
    c   = double(subs(  c,y,z));
end
