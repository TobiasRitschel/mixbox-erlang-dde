function [f, dfx, dfz, d2fx, d2fz, d2fxz] = ...
    logistic_right_hand_side_function(~, x, z, p)
% Parameters
kappa = p.kappa;
K     = p.K;

% Evaluate right-hand side function
f = kappa*x*(1 - z/K);

if(nargout > 1) % Compute Jacobians
    % Jacobian wrt. current states
    dfx = kappa*(1 - z/K);
    
    % Jacobian wrt. delayed states
    dfz = -kappa*x/K;
    
    if(nargout > 2) % Compute Hessians
        % Hessian wrt. states
        d2fx = 0;

        % Hessian wrt. delayed states
        d2fz = 0;

        % Hessian wrt. states and delayed states
        d2fxz = -kappa/K;
    end
end