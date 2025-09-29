function x = compute_steady_state_ddde(x0, f, h, fsolve_opts, varargin)

% Solve residual equation for steady state
x = fsolve(@residual_equation_ddde, x0, fsolve_opts, ...
    f, h, varargin{:});
end

function [res, J] = residual_equation_ddde(x, f, h, varargin)
if(nargout > 1)
    % Evaluate the delayed quantity
    [r, drx] = h(x, varargin{:});

    % Evaluate the delayed state
    z = r;

    % Evaluate the residual function
    [res, dfx, dfz] = f(x, z, varargin{:});

    % Jacobian of the delayed state
    dzx = drx;

    % Jacobian
    J = dfx + dfz*dzx;
else
    % Evaluate the delayed quantity
    r = h(x, varargin{:});

    % Evaluate the delayed state
    z = r;

    % Evaluate the residual function
    res = f(x, z, varargin{:});
end
end