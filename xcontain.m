function [b,gbar] = xcontain(p,x,g,xbar)

problem.solver  = 'fmincon';
problem.options = optimoptions('fmincon', 'Display','none');

problem.objective = @(b) -double(subs(p,x,b*xbar));
problem.nonlcon   = @(b) nonlcon(b,p,x,g,xbar);
problem.lb = 0;
problem.x0 = 0;

[b,fval,flag] = fmincon(problem);

if flag > 0
    gbar = -fval;
else
    b = [];
    gbar = 0;
end

end

function [c,ceq] = nonlcon(b,p,x,g,x1)
    c = double(subs(p,x,b*x1)-g);
    ceq = [];
end
