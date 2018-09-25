function varargout = size(sp,varargin)
% Size of spline functions.

Z = zeros(sp.matdim);

varargout = cell(1,nargout);

[varargout{:}] = size(Z,varargin{:});

end