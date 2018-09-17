function [gbnds,sout1,sout2,info] = upcontain(pi,pb,z,zi,opts)
% Maximizes g subject to the set union containment constraint
%
%   {x: pb(x) <= g} is subset of the union {x: pi(x) >= 0}
%
% where pi is a k-by-1 vector of polynomials.
%
% The set union containment problem can be encoded as multiple piecewise
% intersection problems
%
%   {x: pb(x) <= g} intersecting {pi*(x) <= 0} is subset of {x: p1(x) >= 0}
%   ...
%   {x: pb(x) <= g} intersecting {pi*(x) <= 0} is subset of {x: pk(x) >= 0}
%
% where pi* is the (k-1)-by-1 vector of polynomials with pi*(x) <= 0 if and
% only if pij(x) <= 0 f.a. 1 <= j <= k and j != i.
%
% UPCONTAIN contains PCONTAIN for unions.
%
%% Usage & description
%
% Inputs:
%       -pi:    polynomials to describe united superset
%       -pb:    polynomial to describe subset
%       -z,zi:  Nz-by-1 column vectors of monomials used to specify the 
%               SOS decision variables s(x), si(x) in the Gram matrix form, 
%               s(x)=z(x)'*C*z(x); si(x)=zi(x)'*C*zi(x).
%       -opts:  options for optimization; see GSOSOPTIONS.
%
% Outputs:
%       -gbnds: 1-by-2 vector [glb,gub] providing lower and upper bound on
%               the maximum value of g; will be empty if no feasible lower
%               bound is found.
%       -s0:    k-by-1 vector of multiplier function proving set 
%               containment of pb with glb
%       -sj:    k-by-k matrix of multiplier functions with empty diagonal 
%               proving union containment with glb;
%               s0, sj will be empty if no feasible lower bound is found.
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
x = unique([pi.varname; pb.varname]);

% default options
if ~exist('opts','var') || isempty(opts)
    opts = gsosoptions;
end
% default multiplier z
% See PCONTAIN
if isempty(z)
    % This is roughly based on guidelines from Appendix A.1.1 
    % of the Ph.D. Thesis by Weehong Tan
    zmind = ceil(pi.mindeg/2);
    zmaxd = ceil((pi.maxdeg-pb.maxdeg) / 2);    
    if zmaxd<zmind
        z=1;
    else
        z = monomials(x, zmind:zmaxd);
    end
end
% default multiplier zi
if ~exist('zi','var') || isempty(zi)
    zi = z;
elseif iscell(zi)
    zi = zi{end};
end

% ensure H is column vector
if isrow(pi)
    pi = pi';
else
    assert(iscolumn(pi), 'H must be row or column of polynomials');
end

% number of constraints Hi
k = length(pi);


%% Call GSOSOPT
% ATTENTION! observe positive definition of pi in {x: pi(x) >= 0}.
%
% to solve:
%       max g such that
%       s0, sij in SOS
%       -( -pi + (g-pb)*s0 - pij*sij ) in SOS f.a. i and j != i
%
% GSOSOPT solves minimization of t := -g
%       min t such that
%       s0, sij in SOS
%       pb*s0 + pi + t*s0 + pij*sij in SOS f.a. i and j != i
t   = pvar('g');
s0  = polynomial(zeros(size(pi)));
sj  = polynomial(zeros(length(pi)));
for i=1:k
    s0(i) = sosdecvar(['c' num2str(i)],z);
    for j=1:k
        if i == j, continue;    end
        %else
        sj(i,j) = sosdecvar(['c' num2str(i) num2str(j)],zi);
    end
end
   
sosc = [
    s0 >= 0
    nonzero(sj) >= 0
    (pb*s0 + pi) + t*s0 + sj*pi >= 0
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
    sout2 = subs(sj, dopt);
else
    gbnds = [];
    sout1 = polynomial;
    sout2 = polynomial;
end

end

function a = nonzero(a)
    a = a(:);
    iz = isequal(a,0);
    a = a(~iz);
end
