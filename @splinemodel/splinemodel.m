classdef splinemodel
% Graph-representation of a spline model.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-09-09
% * Changed:    2018-09-09
%
%%    

properties (Access=protected)
    % adjacence matrix
    adj;
    
    % boundary condition matrix
    H;
end
   
properties (Access=public)
    % spline functions
    f;
end
   
properties (SetAccess=private, Dependent)
    % dimension of spline functions
    matdim;
    
    % maximum polynomial degree
    maxdeg;
    
    % minimum polynomial degree
    mindeg;
end

methods (Access=protected)
    sp = combine(sp1,sp2,op);
end
    
methods
    function obj = splinemodel(varargin)
        % Creates a new spline model.
        
        if length(varargin) == 1 && isnumeric(varargin{1})
            f = cell(1,varargin{1});
        else
            f = varargin;
        end
        
        obj.adj = false(length(f));
        obj.H = mpvar('h', size(obj.adj));
        
        obj.H(:,:) = 0;
        obj.f = f;        
    end
        
    function obj = set_adjacent(obj, i, j, h)
        % Sets nodes i and j to be adjacent with boundary condition h s.t.
        %
        %   x in Pi - Pj => h(x) <= 0
        %   x in Pj - Pi => h(x) >= 0
        
        obj.adj(i,j) = true;
        obj.adj(j,i) = true;
        
        obj.H(i,j) = h;
        obj.H(j,i) = -h;
    end
    
    function [b, h] = is_adjacent(obj, i, j)
        % Checks whether nodes i and j are adjacent and returns boundary
        % condition h with
        %
        %   x in Pi - Pj => h(x) <= 0
        %   x in Pj - Pi => h(x) >= 0
        %
        % If the nodes are not adjacent, the boundary condition is zero.
        
        b = obj.adj(i,j);
        h = obj.H(i,j);
    end
    
    function [J, H] = adjacent(obj, i)
        % Returns adjacents j of node i and the respective boundary
        % conditions hij.
        
        idx = 1:count(obj);
        
        J = idx(obj.adj(i,:));
        H = obj.H(i,J);
    end
    
    function [J, H] = paths(obj, i, I)
        % Returns all adjacents j of node i that are in I and their
        % respective boundary conditions hij.
        
        [J, H] = obj.adjacent(i);
        [J, b] = intersect(J, I);
        H = H(b);
    end
    
    function J = adjacents(obj, I)
        % Returns all adjacents j of I that are not already in I.
        
        idx = 1:count(obj);
        
        J = obj.adj(I,:);
           
        J = setdiff(idx(any(J,1)), I);
    end
    
    function H = getH(obj)
        % Returns matrix of boundaries as N-by-N polynomial matrix. No
        % further properties guaranteed.
        
        H = obj.H;
    end
    
    function obj = set.f(obj,value)
        % Set spline functions.
        assert(iscell(value) && length(value) == count(obj), 'Spline functions must be cell array of matching length.');
        
        obj.f = value;
    end
    
    function md = get.matdim(obj)
        md = size(obj.f{end});
    end
    
    function md = get.maxdeg(obj)
        md = max(cellfun(@(c) polynomial(c).maxdeg, obj.f));
    end
    
    function md = get.mindeg(obj)
        md = min(cellfun(@(c) polynomial(c).mindeg, obj.f));
    end
    
    function obj = subsasgn(obj,s,varargin)
        % See SUBSASGN.
        
        if length(s) == 1 && strcmp(s(1).type, '()') && length(s.subs) == 2 && isa(varargin{1}, 'polynomial')
            warning('Use of SP(i,j) for adjacents is deprecated. Use SP.s(i,j) instead.');
            obj = set_adjacent(obj, s.subs{1}, s.subs{2}, varargin{:});
        elseif length(s) == 1 && strcmp(s(1).type, '()')
            obj = subsasgn_idx(obj,s, varargin{:});
        elseif length(s) == 2 && strcmp(s(1).type, '.') && strcmp(s(1).subs, 's') && strcmp(s(2).type, '()')
            obj = set_adjacent(obj, s(2).subs{1}, s(2).subs{2}, varargin{:});
        elseif length(s) == 2 && strcmp(s(1).type, '.') && strcmp(s(1).subs, '->') && strcmp(s(2).type, '()')
            obj = set_adjacent(obj, s(2).subs{1}, s(2).subs{2}, varargin{:});
        elseif length(s) == 2 && strcmp(s(1).type, '.') && strcmp(s(1).subs, '<-') && strcmp(s(2).type, '()')
            error('Use of <- is forbidden for spline systems.')
%             obj = set_adjacent(obj, s(2).subs{1}, s(2).subs{2}, -varargin{:});
        elseif length(s) == 2 && strcmp(s(1).type, '.') && strcmp(s(1).subs, 'f') && strcmp(s(2).type, '()')
            obj = builtin('subsasgn', obj, s, varargin);
        else
            obj = builtin('subsasgn', obj, s, varargin{:});
        end
    end
    
    function out = subsref(obj,s)
        % See SUBSREF.
        
        if length(s) == 1 && strcmp(s.type, '()')
            out = subsref_idx(obj,s);
        else
            out = builtin('subsref', obj, s);
        end
    end
end
    
end

