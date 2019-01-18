function [gbnds,sout1,sout2,info] = spcontain(pa,p2,H,z,zi,opts)
% Maximizes g subject to the spline set containment constraint
%
%   {x:p2(x) <= g} intersecting {x:H(x) <= 0} is subset of {x:pa(x) <= 0}, 
%
% where H is a k-by-1 vector of polynomials with H(x) <= 0 if and only if 
% Hi(x) <= 0 f.a. 1 <= i <= k.
%
% SPCONTAIN extends PWPCONTAIN for multiple boundaries.
%
%% Usage & description
%
% Inputs:
%       -pa:    polynomial to describe piecewise superset
%       -p2:    polynomial to describe subset
%       -H:     polynomials to describe boundary conditions
%       -z,zi:  Nz-by-1 column vectors of monomials used to specify the 
%               SOS decision variables s(x), si(x) in the Gram matrix form, 
%               s(x)=z(x)'*C*z(x); si(x)=zi(x)'*C*zi(x).
%       -opts:  options for optimization; see GSOSOPTIONS.
%
% Outputs:
%       -gbnds: 1-by-2 vector [glb,gub] providing lower and upper bound on
%               the maximum value of g; will be empty if no feasible lower
%               bound is found.
%       -s0:    Multiplier function proving set containment of p2 with glb
%       -si:    k-by-1 vector of multiplier functions proving spline set
%               containment with glb;
%               s0, si will be empty if no feasible lower bound is found.
%       -info:  information structure returned by GSOSOPT
%
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-09-10
% * Changed:    2018-09-10
%
%% See also
%
% See PWPCONTAIN, PCONTAIN
%%

% polynomial variables
x = unique([pa.varname; p2.varname; H.varname]);

% default options
if ~exist('opts','var') || isempty(opts)
    opts = gsosoptions;
end

% ensure H is column vector
if isrow(H)
    H = H';
else
    assert(iscolumn(H), 'H must be row or column of polynomials');
end

% number of constraints Hi
if isempty(H)
    % fall back to polynomial containment problem
    % See PCONTAIN.
    k  = 0;
else
    k = length(H);
end

% origin in cell
if iscell(zi) && all(double(subs(H, x, zeros(size(x)))) <= 0)
    zi = zi{1};
elseif iscell(zi)
    zi = zi{end};
end

%% Call GSOSOPT
% to solve:
%       max g such that
%       s, si in SOS
%       -( pa + (g-p2)*s - H'*si ) in SOS
%
% GSOSOPT solves minimization of t := -g
%       min t such that
%       s, si in SOS
%       p2*s - pa + t*s + H'*si in SOS
t   = pvar('g');
s0  = sosdecvar('c',z);
si  = sosmdecvar('d',zi,k);


sosc = [
    s0 >= 0
    si >= 0
    (p2*s0 - pa) + t*s0 + nonempty(H'*si) >= 0
];

gopts = opts;
gopts.minobj = -opts.maxobj;
gopts.maxobj = -opts.minobj;
if strcmp(opts.display,'on')
    gopts.display = 'pcontain'; % Undocumented syntax for proper display
end

[info,dopt] = gsosopt(sosc,x,t,gopts);

if ~isempty(info.tbnds)
    gbnds = -info.tbnds([2 1]);
    
    sout1 = subs(s0, dopt);
    sout2 = subs(si, dopt);
else
    gbnds = [];
    sout1 = polynomial;
    sout2 = polynomial;
end

end

function a = nonempty(a)
    if isempty(a)
        a = 0;
    end
end