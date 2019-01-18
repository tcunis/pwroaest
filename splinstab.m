function varargout = splinstab(f, H, x, Q, Adj, opts)
% Performs a linear stability analysis for a switched polynomial system,
%
%           |
%           |
%   fi(x)   |   fj(x)
%           |
%           |
%       hij(x) < 0
%
% for all 1 <= i < j <= k, about the equilibrium point x = 0.
%
%% Usage & description
%
%   V = pwlinstab(f,H,x)
%   [V1,...,Vk] = pwlinstab(f,H,x)
%   [...] = pwlinstab(...,Q)
%   [...] = pwlinstab(...,Q,opts)
%
% Inputs:
%       -f:   k-by-1 cell of polynomials vector fields
%       -H:   k-by-1 cell of boundary conditions (scalar fields);
%             f,H correspond to the k spline domains with x=0 at the 
%             boundary.
%       -x:   state-space vector as PVAR
%       -Q:   postive definite parameter in the Lyapunov equations; must be
%             either scalar or square matrix of the size of |x|, or empty. 
%             [default = 1e-6]
%       -opts:  SOS options structure; see SOSOPTIONS.
%
% Outputs:
%       -V:   common quadratic Lyapunov function, V = x'*P*x
%       -V1,V2   multiple quadratic Lyapunov functions, Vi = x'*Pi*xi; 
%                with hij(x) = 0 -> Vi(x) = Vj(x) for all 1 <= i < j <= k.
%
% If the linearization is not stable then V and P are returned as empty.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-09-20
% * Changed:    2019-01-11
%
%%

assert(nargout == 1 || nargout == length(f), 'Undefined number of outputs (%g).', nargout);


if ~exist('Q', 'var') || isempty(Q)
    Q = 1e-6;
end
if ~exist('Adj', 'var') || isempty(Adj)
    Adj = [];
end
if ~exist('opts', 'var')
    opts = sosoptions;
end

k = length(f);

varargout = cell(1,nargout);

A  = cell(k,1);
ev = cell(k,1);

% linearize: xdot = Ai*x
for i=I
    A{i}  = plinearize(f{i},x);
    ev{i} = eig(A{i});
end

if max(real(vertcat(ev{:}))) >= 0
    % if any Ai is unstable
    % nothing to do
elseif nargout == 1
    % if A1, A2 are stable solve LMI
    %   A1'*P + P*A1 + S1 < 0
    %   A2'*P + P*A2 + S2 < 0
    % for common Lyapunov function
    [V,P] = sosdecvar('p', x);
    
    sosc = polynomial(zeros(k+1,1)) == 0;
    
    for i=1:k
        sosc(i) = x'*(A{i}'*P + P*A{i})*x - sum(H{i}) <= -x'*Q*x;
    end

    sosc(end) = V >= x'*Q*x;

    [info,dopt] = sosopt(sosc, x, opts);
    
    if info.feas
        P = subs(P, dopt);
        varargout{1} = x'*P*x;
    end
else
    % if A1, A2 are stable solve LMI
    %   A1'*P + P1*A1 < 0
    %   A2'*P + P2*A2 < 0
    % for multiple Lyapunov function
    V = cell(k,1);
    P = cell(k,1);
    
    soscP = polynomial(zeros(k,1)) == 0;
    soscV = polynomial(zeros(k,1)) == 0;
    
    % continuity decision variables
    z  = monomials(x, 0:2);
    
    for i=1:k
        [V{i},P{i}] = sosdecvar(['p' num2str(i)], x);
        
        soscP(i) = x'*(A{i}'*P{i} + P{i}*A{i})*x <= x'*Q*x;
        soscV(i) = V{i} >= x'*Q*x;
    end
     
    soscH = sproav_continuity(V, H, Adj, z);
    
    sosc = [soscP; soscV; soscH];
    
    [info,dopt] = sosopt(sosc, x, opts);
    
    if info.feas
        for i=1:k
            varargout{i} = subs(V{i},dopt);
        end
    end
end

end

    