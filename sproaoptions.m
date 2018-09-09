classdef sproaoptions < pwroaoptions
% Options for spline ROA estimation.
%
%% Usage & description
%
%   opts = sproaoptions(SP, x, ...)
%
% with inputs
%       -SP:  spline-model
%       -x:   state-space vector as PVAR
%
%   Name, Value:
%       -xi: variables of boundary condition as PVAR; must be subset of the
%            state-space variables. [default =  opts.x]
%       -zi: Monomials for boundary multiplier in gamma-s2-si step of V-s 
%            iteration. Specifically, zi is a column vector of monomials 
%            used to specify si(xi) in the Gram matrix form, 
%            si(xi)=zi(xi)'*C*zi(xi). [default =  monomials(xi, 1:2) ]
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
% See PWROAOPTIONS, ROAOPTIONS
%%

properties
    %x -- inherited from ROAOPTIONS
    SP;
    %zi -- inherited from PWROAOPTONS
    %zV -- inherited from ROAOPTIONS
    %zVi -- inherited from PWROAOPTIONS
end

methods
    function opt = sproaoptions(SP, x, varargin)
        opt@pwroaoptions(SP.f{1}, SP.f(2:end), SP.getH, x, varargin{:});
        
        opt.SP = SP;
    end
end

end