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
% * Changed:    2019-01-11
%
%% See also
%
% See PWROAEST, SPROAOPTIONS, PCONTAIN, PWPCONTAIN
%%

% information from options
p  = roaopts.p;
zV = roaopts.zVi;
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
debug = roaopts.debug;
Vin = roaopts.Vi0;
Ik = roaopts.Ik;

% Vdeg = zV.maxdeg;
Nsteps = NstepBis;

% initialize storage
c0 = cell(Nsteps,1);
iter= struct('V',c0,'beta',c0,'gamma',c0,'Ik',c0,'s0',c0,'s',c0,'si',c0,'it2',c0,'count',c0,'time',c0);

% Lyapunov functions
V = cell(size(zV));
if length(zV) > 1
    c = @(i) i;
else
    c = @(~) 1;
end

if isempty(Ik)
    % find active domain at origin
    [~,Ik] = double(subs(SP, x, zeros(size(x))));
end

%% Run V-s iteration
fprintf('\n---------------Beginning spline V-s iteration\n');
biscount = 0;
for k1=1:NstepBis
    time = [];
    time_tot = tic;
    
    time_vstep = tic;
    %======================================================================
    % Find V step:
    % Hold s1, s2, si, b, g fixed and solve the Gamma Step and Beta Step
    % for Appropiate Lyapunov function with an additional constraint
    % V-L1 in SOS
    %======================================================================
    if k1==1
        if ~isempty(Vin)
            V(c(Ik)) = Vin;
        elseif length(Ik) == 1
            % construct Lyap function from linearization of active domain
            V{c(Ik)} = linstab(SP.f{Ik},x,Q);
        else
            H = cell(count(SP),1);
            A = false(count(SP));
            for i=Ik, [J,H{i}] = paths(SP,i,Ik); A(i,Ik) = ismember(Ik,J); end
            [V{c(Ik)}] = splinstab(SP.f(Ik),H(Ik),x,Q,A(Ik,Ik));
        end
        
    elseif length(Ik) == 1
        % local V-s problem
        V{c(Ik)} = roavstep(SP.f{Ik},p,x,zV,b,g,sb{c(Ik)},s{Ik},L1,L2,sopts);
        if isempty(V{c(Ik)})
            if strcmp(display,'on')
                fprintf('local V-step infeasible at iteration = %d\n',k1);
            end
            break;
        end
    else
        [V{c(Ik)}] = sproavstep(SP.f(Ik),H(Ik),p,x,zV,beta(c(Ik)),g,sb(c(Ik)),s(Ik,:),A(Ik,Ik),L1,L2,roaopts);
        if isempty(V{c(1)})
            if strcmp(display,'on') && length(V) == 1
                fprintf('common V-step infeasible at iteration = %d\n',k1);
            elseif strcmp(display,'on')
                fprintf('multiple V-step infeasible at iteration %d\n',k1);
            end
            break;
        end
    end
    time.vstep = toc(time_vstep);
    
    % boundaries, adjacents, & multipliers
    s = cell(count(SP),2);
    H = cell(count(SP),1);
    A = false(count(SP));
    
    time.gpre = zeros(1,count(SP)-length(Ik)+1);
    time.gmin = zeros(1,count(SP)-length(Ik)+1);
    time.next = zeros(1,count(SP)-length(Ik)+1);
    
    quota = 0;

    for k2=1:count(SP)-length(Ik)+1
        time_gpre = tic;
        %==================================================================
        % Pre Gamma Step: Solve the problem max g s.t. for all i in Ik
        % {x:V(x) <= gamma} minus {x:hij(x) >= 0} f.a. j in J
        % is contained in {x:grad(V)*fi < 0}
        % with J = Ik - {i} = paths[k](i)
        %==================================================================
        gopts.maxobj = gammamax;
        gpre = zeros(size(Ik));
        for i=Ik
            if k2 > 1 && (i ~= next) && ~is_adjacent(SP,i,next)
                gpre(Ik==i) = gpre2(Ik(2:end)==i);
                continue;
            end
            
            % else:
            [J,H{i}] = paths(SP,i,Ik);
            [gbnds,s{i,:}] = spcontain(jacobian(V{c(i)},x)*SP.f{i}+L2,V{c(i)},H{i},z2,zi,gopts);
            if isempty(gbnds)
                if strcmp(display,'on')
                    fprintf('pre gamma step for domain %d infeasible at iteration = %d-%d.\n', i, k1, k2);
                end
                break;
            end
            gpre(Ik==i) = gbnds(1);
            
            if gbnds(2) < gopts.maxobj
                gopts.maxobj = gbnds(2);
            end

            % adjacents matrix for multiple V-step
            A(i,Ik) = ismember(Ik,J);

            quota = quota + length(H{i}) + 1;
            
            if strcmp(debug,'on')
                fprintf('debug: i = %d \t gpre = %4.6f\n', i, gpre(Ik==i));
            end
        end
        time.gpre(k2) = toc(time_gpre);
        % sort & store
        [gpre2,is] = sort(gpre);
        Ik = Ik(is);
        % take minimum
        gpre = min(gpre);
        
        time_gmin = tic;
        %==================================================================
        % Min Gamma Step: Solve the problem max g s.t. for all i in adj(Ik)
        % {x:V(x) <= gamma} is contained in the union {x:hij(x) >= 0} f.a.
        % j in J
        % with J = paths[k](i)
        %==================================================================
        I = adjacents(SP, Ik);
%         gopts.maxobj = gammamax;
        opts = optimset('Display','notify');
        dist = zeros(size(I));
        for i=I
%             if k2 > 1 && any(I2==i)
%                 dist(I==i) = dist2(I2==i);
%                 continue;
%             end
            
            % else:
            [J,Hj] = paths(SP,i,Ik);
            [~,Hi] = paths(SP,i,I);
            
            distV = zeros(c(length(J)),1);
            for j=c(J)
                [gbnds,~] = spdistance(Hj,V{j},Hi,opts);
                if isempty(gbnds)
                    if strcmp(display,'on')
                        fprintf('distance step for domain %d infeasible at iteration = %d-%d.\n', i, k1, k2);
                    end
                    break;
                end
                distV(c(J)==j) = gbnds(1);
            end
            dist(I==i) = min(distV);
            
            if strcmp(debug,'on')
                fprintf('debug: i = %d \t dist = %4.6f\n', i, dist(I==i));
            end
        end
        time.gmin(k2) = toc(time_gmin);
        % store
%         I2 = I;
%         dist2 = dist;
        % find minimum
        [gmin,ig] = min(dist);
        
        % print results and proceed
        if strcmp(display,'on') && length(Ik) < length(SP)
            fprintf('iteration = %d-%d\t gpre = %4.6f\t gmin = %4.6f\t Ik = (%s)\t J = (%s)\n',k1,k2,gpre,gmin,num2str(Ik),num2str(I));
        end
        
        time_next = tic;
        
        if gpre > gmin
            next = I(ig);                
            Ik = [next Ik];
            
            if length(V) > 1
                % find viable Vnext
                [I,H{next}] = paths(SP,next,Ik);
                A(next,Ik) = ismember(Ik,I);
                % symmetric adjacents
                for i=I
                    [J,H{i}] = paths(SP,i,Ik);
                    A(i,Ik) = ismember(Ik,J);
                end
                
                % use quadratic Lyapunov-candidate in first iteration
                z = zV{next};
                V{next} = sproavnext(V(Ik),H(Ik),x,A(Ik,Ik),z,zi);
                if isempty(V{next})
                    if strcmp(display,'on')
                        fprintf('next V-step infeasible at iteration %d-%d\n',k1,k2);
                    end
                    break;
                end
            end
        else
            g = gpre;
            break;
        end
        
        time.next(k2) = toc(time_next);
    end
    
    time_beta = tic;
    
    % pseudo-radii & multipliers
    beta = zeros(length(V),1);
    sb = cell(length(V),1);
    
    for i=c(Ik)
        %======================================================================
        % Beta Step: Solve the problem max b s.t.
        % {x:p(x) <= b} is contained in {x:V(x)<=g}
        %======================================================================
        gopts.maxobj = betamax;
        [bbnds,sb{i}] = pcontain(V{i}-g,p,z1,gopts);
        if isempty(bbnds)
            if strcmp(display,'on')
                fprintf('beta step infeasible at iteration = %d\n', k1);
            end
            break;
        end
        beta(i) = bbnds(1);
    end
    
    time.beta = toc(time_beta);
    
    b = min(beta(c(Ik)));

    % print results and store iteration data
    if strcmp(display,'on')
        fprintf('iteration = %d\t beta = %4.6f\t gamma = %4.6f\t Ik = (%s)\n',k1,b,g,num2str(Ik));
    end
    time.gpre = sum(time.gpre);
    time.gmin = sum(time.gmin);
    time.next = sum(time.next);
    time.tot  = toc(time_tot);
    iter(k1).V      = V;
    iter(k1).beta   = b;
    iter(k1).gamma  = [g gpre gmin];
    iter(k1).Ik     = Ik;
    iter(k1).s0     = sb;
    iter(k1).s      = s(Ik,1);
    iter(k1).si     = s(Ik,2);
    iter(k1).it2    = k2;
    iter(k1).time   = time;
    iter(k1).count  = quota;
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