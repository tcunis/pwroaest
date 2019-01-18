function varargout = spcontour(p,v,H,domain,linespec,npts,var)
% Plots contours of p(x,y) within domain specified by Hi(x,y) <= 0 f.a. i.
%
%% About
%
% * Author:     Torbjoern Cunis, PJS
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2019-01-11
% * Changed:    2019-01-11
%
%% Based on
%
% See PCONTOUR
%%

if nargin==3    
    domain=[];
    linespec=[];
    npts=[];
    var=[];
elseif nargin==4
    linespec=[];
    npts=[];
    var=[];
elseif nargin==5
    npts=[];
    var=[];
elseif nargin==6
    var=[];
end

% Default contour
if isempty(v)
    v=1;
end
lv = length(v);

% Default linespec
if isempty(linespec)
    linespec = 'b';
end

% Default npts
if isempty(npts)    
    Nx = 100;
    Ny = 100;
else
    Nx = npts(1);
    if length(npts)==1
        Ny = Nx;
    else
        Ny = npts(2);
    end
end

% Define variables as chars
if isempty(var)
    x = p.varname{1};
    y = p.varname{2};
elseif ispvar(var) || iscellstr(var)
    if ispvar(var)
        var = char(var);
    end
    x = var{1};
    y = var{2};
else
    error('var must be a vector of pvars');
end

if isempty(domain)
    domain = [-1 1 -1 1];
end    
    
% Plot contour
xg = linspace(domain(1),domain(2),Nx);
yg = linspace(domain(3),domain(4),Ny);
[xg,yg] = meshgrid(xg,yg);
pgrid = double(subs(p,{x; y},[xg(:)'; yg(:)']));
Hgrid = double(subs(H,{x; y},[xg(:)'; yg(:)']));

pgrid(any(Hgrid>0,1)) = NaN;
pgrid = reshape(pgrid,size(xg));

if lv==1
    % Single contour syntax for contour function
    v = [v v];
end

if nargout==0
    contour(xg,yg,pgrid,v,linespec);    
else
    [C,h]=contour(xg,yg,pgrid,v,linespec);     
    varargout = {C,h};
end
xlabel(x)
ylabel(y)

