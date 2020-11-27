function varargout = pcontour(f,obj,varargin)
% Plots contours of function family over spline domains.
%
% If no function is given, the spline domain boundaries are plotted.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:tcunis@umich.edu>
% * Created:    2019-07-09
% * Changed:    2019-07-09
%
%% See
%
% See PCONTOUR, SPCONTOUR
%%

isholdon = ishold;
hold on

if isa(f,'splinemodel')
    % no function given
    varargin = [{obj} varargin];
    obj = f;
    f = [];
end

if isempty(f)
    % Plot domain boundaries
    H = triu(obj.getH);
    H = H(~isequal(H,0));
    
    x = polynomial(H.varname);
    
%     C = cell(1,length(H));
%     h = cell(1,length(H));
    
    for i=1:length(H)
%         [C{i},h{i}]=
        pcontour(H(i)+1e-12*(x'*x), varargin{:});
    end
    
%     h = [h{:}];
else
    % Plot functions over domains
    idx = 1:count(obj);
    I = idx(cellfun(@(p) ~isempty(p), f));
    
%     C = cell(1,length(idx));
%     h = cell(1,length(idx));
    
    for i=I
        [~,Hi] = paths(obj,i,I);
        
%         [C{i},h{i}]=
        spcontour(f{i},varargin{1},Hi,varargin{2:end});
    end
end

if ~isholdon
    hold off
end

% if nargout > 0
%     varargout = {C,h};
% end
varargout = {};

end