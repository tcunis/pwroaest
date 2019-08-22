function varargout = sproavstep(f,H,p,x,z,beta,gamma,sb,s,Adj,L1,L2,roaopts,Hnew)
% Solves the V-s step of the spline ROA iteration.
%
%% Usage & description
%
%   [V,c] = sproavstep(f,H,p,x,z,beta,gamma,s0,s)
%   [...] = sproavstep(...,L1,L2,opts)
%
% Solves for a common Lyapunov function V which satisfies the SOS
% constraints for the previously covered spline domains while holding all 
% other variables fixed.
%
% Inputs:
%       -f:   k-by-1 cell of polynomials vector fields
%       -H:   k-by-1 cell of boundary conditions (scalar fields);
%             f,H correspond to the k previously covered spline domains.
%       -p:   shape function (scalar field)
%       -x:   state-space vector as PVAR
%       -z:   Nz-by-1 column vector of monomials; specifies the Lyapunov
%             function decision variable in the vector form V(x) = c'*z(x).
%       -beta:  level set of shape function
%       -gamma: level set of Lyapunov function
%       -s0:  multiplier for beta-step 
%       -s:   k-by-2 cell of multipliers for gamma-step (first column: 
%             Lyapunov gradient; second column: boundary condition)
%       -L1:  epsilon 1 (double or scalar field); enforces strict positive 
%             definiteness of the Lyapunov function.
%       -L2:  epsilon 2 (double or scalar field); enforces strict negative
%             definiteness of the Lyapunov gradients.
%       -opts:  SOS options structre; see SOSOPTIONS.
%
% Outputs:
%       -V: Lyapunov function satisfying the SOS constraints; will be empty
%           if no feasible solution is found.
%       -c: Nz-by-1 vector of coefficients of the returned Lyapunov
%           function.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-09-10
% * Changed:    2019-01-11
%
%% See also
%
% See SPROAEST, PWROAVSTEP
%%

if isa(roaopts, 'sosoptions')
    opts = roaopts;
    
    zi = 1;
else
    opts = roaopts.sosopts;
    
    zi = roaopts.zi{end};
end

if ~exist('Hnew','var')
    Hnew = H;
end

% number of domains covered
k = length(f);

assert(length(H) == k && size(s,1) == k, ...
       'Number of functions, boundaries, and multipliers must be equal to number of domains.' ...
);

if length(z) == 1
    % common Lyapunov function
    varargout = cell(1,2);

    % Lyapunov decision variable
    [V,c] = polydecvar('c',z{:});

    %% Common V-s feasibility problem
    % prepare constraint variable
    sosconstr = polynomial(zeros(k+2,1)) == 0;

    % V-L1 in SOS
    sosconstr(1) = V >= L1;

    % {x: p(x) <= b} is contained in {x: V(x) <= g}
    sosconstr(2) = -((V-gamma) + sb{:}*(beta-p)) >= 0;

    % {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
    gradV = jacobian(V,x);

    for i=1:k
        % -( pa + (g-p2)*s - H'*si ) in SOS
        sosconstr(2+i) = -(gradV*f{i} + L2 + s{i,1}*(gamma-V) - H{i}'*s{i,2}) >= 0;
    end

    % solve problem
    [info,dopt] = sosopt(sosconstr,x,opts);

    % output
    if info.feas
        varargout{1}   = subs(V,dopt);
        varargout{end} = subs(c,dopt);
    end
else
    % multiple Lyapunov functions
    varargout = cell(1,2*k);
    
    % Lyapunov decision variable
    V = cell(k,1); c = cell(k,1);
    
    % prepare constraint variable
    soscV = polynomial(zeros(k,1)) == 0;
    soscB = cell(k,1); %polynomial(zeros(k,1)) == 0;
    soscG = cell(k,1); %polynomial(zeros(k,1)) == 0;
    
    %% Multiple V-s feasibility problem
    for i=1:k
        [V{i},c{i}] = polydecvar(['c' num2str(i)],z{i});
        
        % Vi-L1 in SOS
        soscV(i) = V{i} >= L1;
        
        if isempty(sb{i})
            % compute Lyapunov-function for new domain
            continue
        end
        
        % {x: p(x) <= bi} is contained in {x: Vi(x) <= g}
        soscB{i} = -((V{i}-gamma) + sb{i}*(beta(i)-p)) >= 0;

        % {x: Vi(x) <= g} is contained in {x: grad(Vi)*fi < 0}
        gradVi = jacobian(V{i},x);

        % -( pa + (g-p2)*s - H'*si ) in SOS
        soscG{i} = -(gradVi*f{i} + L2 + s{i,1}*(gamma-V{i}) - H{i}'*s{i,2}) >= 0;
    end
    
    % continuity
    soscH = sproav_continuity(V, Hnew, Adj, zi);
    
    sosconstr = vertcat(soscV, soscB{:}, soscG{:}, soscH);
    
    % solve problem
    [info,dopt] = sosopt(sosconstr,x,opts);
    
    % output
    if info.feas
        for i=1:k
            varargout{i}   = subs(V{i},dopt);
            varargout{k+i} = subs(c{i},dopt);
        end
    end
    
end
