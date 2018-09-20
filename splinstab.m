function varargout = splinstab(f, H, x, Q, opts)
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
% * Changed:    2018-09-20
%
%%

assert(nargout == 1 || nargout == length(f), 'Undefined number of outputs (%g).', nargout);


if ~exist('Q', 'var') || isempty(Q)
    Q = 1e-6;
end
if ~exist('opts', 'var')
    opts = sosoptions;
end

k = length(f);

varargout = cell(1,nargout);

A  = cell(k,1);
ev = cell(k,1);

% linearize: xdot = Ai*x
for i=1:k
    A(i)  = plinearize(f{i},x);
    ev(i) = ev(A{i});
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
    
    sosc = polyconstr(k+1);
    sosc(1) = V >= x'*Q*x;
    
    for i=1:k
        sosc(1+i) = x'*(A{i}'*P + P*A{i})*x - H{i} <= -x'*Q*x;
    end
    
    [info,dopt] = sosopt(sosc, x, opts);
    
    if info.feas
        P = subs(P, dopt);
        varargout{1} = x'*P*x;
    end
else
    error('Multiple Lyapunov functions not yet implemented for splines.');
end

end

    