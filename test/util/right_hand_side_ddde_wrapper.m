function [rhs, J] = right_hand_side_ddde_wrapper ...
    (inp, t, p, varargin)

% Evaluate right-hand side function
rhs = right_hand_side_ddde ...
    (t, inp, p, varargin{:});

if(nargout > 1)
    % Evaluate Jacobian
    J = jacobian_ddde ...
        (t, inp, p, varargin{:});
end