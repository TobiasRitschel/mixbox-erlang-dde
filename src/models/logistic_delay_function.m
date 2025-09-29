function [h, dhx, d2hx] = logistic_delay_function(x, p, varargin) %#ok
% Delayed state
h = x;

if(nargout > 1)
    % Jacobian wrt. states
    dhx = eye(numel(x));

    if(nargout > 3)
        % Hessian wrt. states
        d2hx = 0;
    end
end