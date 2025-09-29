function [f, df, d2f] = right_hand_side_function_wrapper(inp, p, fun)

% Dimensions
nx = p.nx;
nz = p.nz;

% Extract variables
x = inp(      1:nx );
z = inp(nx + (1:nz));

if(nargout > 2)
    % Evaluate right-hand side function
    [f, dfx, dfz, d2fx, d2fz, d2fxz] = ...
        fun([], x, z, p);

    for i = 1:nx
    % Create Hessian
    d2f(:, :, i) = [
        d2fx(:, :, i),  d2fxz(:, :, i);
        d2fxz(:, :, i)', d2fz(:, :, i)]; %#ok
    end
elseif(nargout > 1)
    % Evaluate right-hand side function
    [f, dfx, dfz] = ...
        fun([], x, z, p);
else
    % Evaluate right-hand side function
    f = fun([], x, z, p);
end

if(nargout > 1)
    % Create Jacobian
    df = [dfx, dfz];
end