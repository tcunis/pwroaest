function [x,I] = double(sp)
% Determines active nodes I and returns double value of their
% functions.
%
% See DOUBLE.

idx = 1:count(sp);

I = all(double(sp.H) <= 0, 2);
x = double([sp.f{I}]);

I = idx(I);

end