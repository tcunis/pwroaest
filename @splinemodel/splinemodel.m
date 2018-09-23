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

methods (Access=protected)
    sp = combine(sp1,sp2,op);
end
    
methods
    function obj = splinemodel(varargin)
        % Creates a new spline model.
        
        if length(varargin) == 1 && isnumeric(varargin{1})
            obj.f = cell(1,varargin{1});
        else
            obj.f = varargin;
        end
        
        obj.adj = false(length(obj.f));
        obj.H = mpvar('h', size(obj.adj));
        
        obj.H(:,:) = 0;
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
        
        idx = 1:length(obj);
        
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
        
        idx = 1:length(obj);
        
        J = obj.adj(I,:);
           
        J = setdiff(idx(any(J,1)), I);
    end
    
    function H = getH(obj)
        % Returns matrix of boundaries as N-by-N polynomial matrix. No
        % further properties guaranteed.
        
        H = obj.H;
    end
    
    function obj = subsasgn(obj,s,varargin)
        % See SUBSASGN.
        
        if length(s) == 1 && strcmp(s.type, '()')
            obj = obj.set_adjacent(s.subs{1}, s.subs{2}, varargin{:});
        else
            obj = builtin('subsasgn', obj, s, varargin);
        end
    end
end
    
end

