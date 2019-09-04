function [flag,msg] = assertIsWellFormed(sp)
% Checks that a spline model is well formed according definition.
%
%% Usage & Description
%
%   assertIsWellFormed(sp)
%
% Throws an assertion failed exception if sp is not well formed. 
%
%   [flag,msg] = assertIsWellFormed(sp)
%
% Instead of an exception, returns well-formed flag and message with
%
%   -flag == 0  : Spline model is well formed;
%   -flag == 1  : Spline model has self-jumps;
%   -flag == 2  : Spline model has jump w/ empty condition.
%
%%

[flag,msg,args] = isWellFormed(sp);

if nargout > 0
    % Return flag & message
    msg = sprintf(msg, args{:});
    
else
    % Throw exception if flag > 0
    assert(flag == 0, msg, args{:});
    
end

end


function [flag,msg,args] = isWellFormed(sp)
% Private function.

    idx = 1:count(sp);

    self_jump = diag(sp.adj);
    
    if any(self_jump)
        % self-jump
        flag = 1;
        msg  = 'Domain(s) %s have unpermitted self-jumps.';
        args = {num2str(idx(self_jump))};
        return
    end
    
    empty_jump = isequal(sp.H(sp.adj), 0);
    
    if any(empty_jump)
        % empty switching condition
        Hfr = repmat(idx', 1, count(sp));   Hfr = Hfr(sp.adj);
        Hto = repmat(idx,  count(sp), 1);   Hto = Hto(sp.adj);
        
        flag = 2;
        msg  = 'Jump(s) from [%s] to [%s] have empty conditions.';
        args = {num2str(Hfr(empty_jump)'), num2str(Hto(empty_jump)')};
        return
    end
    
    %else
    flag = 0;
    msg  = '';
    args = {};
    
end
