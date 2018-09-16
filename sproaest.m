function [beta,V,gamma,Ik,iter] = sproaest(SP,x,roaopts)
% Estimates lower bound of spline region of attraction.
%
%% Usage & description
%
%   [beta,V,gamma,Ik,iter] = sproaest(SP, x, roaopts)
%
% Estimates the lower bound of the region of attraction of the spline model
% SP about the equilibrium point x = 0.
%
% Inputs:
%       -SP:  spline model
%       -x:   state-space vector as PVAR
%       -ropts: options for spline ROA estimation; see SPROAOPTIONS.
%
% Outputs:
%       -beta:  maximum beta
%       -V:     Lyapunov function corresponding to beta
%       -gamma: Lyapunov level set corresponding to beta
%       -Ik:    set of spline domains covered by {V <= gamma}
%       -iter:  structure with fields V, beta, gamma, s, and time; 
%               contains the results for each iteration.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-09-09
% * Changed:    2018-09-09
%
%% See also
%
% See PWROAEST, SPROAOPTIONS, PCONTAIN, PWPCONTAIN
%%

% information from options
p  = roaopts.p;
zV = roaopts.zV; %TODO: multiple Lyapunov functions
z1 = roaopts.z1;
z2 = roaopts.z2;
zi = roaopts.zi;
L2 = roaopts.L2;
L1 = roaopts.L1;
Q  = roaopts.Q;
NstepBis = roaopts.NstepBis;
sopts = roaopts.sosopts;
gopts = roaopts.gsosopts;
gopts.minobj = 0;
gammamax = roaopts.gammamax;
betamax = roaopts.betamax;
display = roaopts.display;
Vin = roaopts.Vi0;

% Vdeg = zV.maxdeg;
Nsteps = NstepBis;

% initialize storage
c0 = cell(Nsteps,1);
iter= struct('V',c0,'beta',c0,'gamma',c0,'s0',c0,'s',c0,'si',c0,'time',c0);

% Lyapunov functions
% V = cell(size(zV));

Ik = 1; %TODO: find active domain at origin

%% Run V-s iteration
fprintf('\n---------------Beginning spline V-s iteration\n');
biscount = 0;
for k1=1:NstepBis
    tic;
    
    %======================================================================
    % Find V step:
    % Hold s1, s2, si, b, g fixed and solve the Gamma Step and Beta Step
    % for Appropiate Lyapunov function with an additional constraint
    % V-L1 in SOS
    %======================================================================
    if k1==1
        if ~isempty(Vin)
            V = Vin;
        elseif length(Ik) == 1
            % Construct Lyap function from linearization of active domain
            V=linstab(SP.f{Ik},x,Q);
        else
            V=splinstab();
        end
        
    elseif length(Ik) == 1
        % local V-s problem
        [V,~] = roavstep(SP.f{Ik},p,x,zV,b,g,s0,s{1},L1,L2,sopts);
        if isempty(V)
            if strcmp(display,'on')
                fprintf('local V-step infeasible at iteration = %d\n',k1);
            end
            break;
        end
    else
        [V,~] = sproavstep(SP.f(Ik),H,p,x,zV,b,g,s0,s,L1,L2,sopts);
        if isempty(V)
            if strcmp(display,'on') %&& length(V) == 1
                fprintf('common V-step infeasible at iteration = %d\n',k1);
            end
            break;
        end
    end
    
    for k2=1:length(SP)-length(Ik)+1
        %==================================================================
        % Pre Gamma Step: Solve the problem max g s.t. for all i in Ik
        % {x:V(x) <= gamma} minus {x:hij(x) >= 0} f.a. j in J
        % is contained in {x:grad(V)*fi < 0}
        % with J = Ik - {i} = paths[k](i)
        %==================================================================
        gopts.maxobj = gammamax;
        gpre = zeros(size(Ik));
        s = cell(length(Ik),2);
        H = cell(length(Ik),1);
        for i=Ik
            [~,H{i}] = paths(SP,i,Ik);
            [gbnds,s{Ik==i,:}] = spcontain(jacobian(V,x)*SP.f{i}+L2,V,H{i},z2,zi,gopts);
            if isempty(gbnds)
                if strcmp(display,'on')
                    fprintf('pre gamma step for domain %d infeasible at iteration = %d-%d.\n', i, k1, k2);
                end
                break;
            end
            gpre(Ik==i) = gbnds(1)
        end
        gpre = min(gpre);
            
        %==================================================================
        % Min Gamma Step: Solve the problem max g s.t. for all i in adj(Ik)
        % {x:V(x) <= gamma} is contained in the union {x:hij(x) >= 0} f.a.
        % j in J
        % with J = paths[k](i)
        %==================================================================
        I = adjacents(SP, Ik);
        gopts.maxobj = gammamax;
        dist = zeros(size(I));
        for i=I
            [~,Hi] = paths(SP,i,Ik);
            [gbnds,~] = upcontain(Hi,V,[],zi,gopts);
            if isempty(gbnds)
                if strcmp(display,'on')
                    fprintf('distance step for domain %d infeasible at iteration = %d-%d.\n', i, k1, k2);
                end
                break;
            end
            dist(I==i) = gbnds(1)
        end
        [gmin,ig] = min(dist);
        
        % print results and proceed
        if strcmp(display,'on') && length(Ik) < length(SP)
            fprintf('iteration = %d-%d\t gpre = %4.6f\t gmin = %4.6f\t Ik = (%s)\t J = (%s)\n',k1,k2,gpre,gmin,num2str(Ik),num2str(I));
        end
        if gpre > gmin
            next = I(ig);
            Ik = [Ik next];
        else
            g = gpre;
            break;
        end
    end
    
    %======================================================================
    % Beta Step: Solve the problem max b s.t.
    % {x:p(x) <= b} is contained in {x:V(x)<=g}
    %======================================================================
    gopts.maxobj = betamax;
    [bbnds,s0] = pcontain(V-g,p,z1,gopts);
    if isempty(bbnds)
        if strcmp(display,'on')
            fprintf('beta step infeasible at iteration = %d\n', k1);
        end
        break;
    end
    b = bbnds(1);

    % print results and store iteration data
    if strcmp(display,'on')
        fprintf('iteration = %d\t beta = %4.6f\t gamma = %4.6f\t Ik = (%s)\n',k1,b,g,num2str(Ik));
    end
    iter(k1).V      = V;
    iter(k1).beta   = b;
    iter(k1).gamma  = [g gpre gmin];
    iter(k1).s0     = s0;
    iter(k1).s1     = s{:,1};
    iter(k1).si     = s{:,2};
    iter(k1).time   = toc;
    biscount = biscount+1;
end
if strcmp(display,'on')
    fprintf('---------------Ending V-s iteration.\n');
end
    
%% Outputs
iter(biscount+1:end) = [];
[~, idx] = max([iter.beta]);
beta  = iter(idx).beta;
V     = iter(idx).V;
gamma = iter(idx).gamma(1);

end