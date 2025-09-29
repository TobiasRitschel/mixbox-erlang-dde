function [R, dR] = residual_equations(xkp1, tkp1, xk, rkmjp1, alphatkmsj, dt, p, f, h)
% Compute Jacobian?
compute_jacobian = (nargout > 1);

if(compute_jacobian)
    % Delayed quantities
    [rkp1, drkp1] = h(xkp1, p);
else
    % Delayed quantities
    rkp1 = h(xkp1, p);
end

% Compute delayed states
zkp1 = sum(alphatkmsj.*[rkp1, rkmjp1], 2)*dt;

if(compute_jacobian)
    % Right-hand side function
    [fkp1, dfxkp1, dfzkp1] = f(tkp1, xkp1, zkp1, p);

    % Identity matrix
    I = eye(size(dfxkp1));

    % Sensitivities of delayed states
    dzkp1 = alphatkmsj(:, 1).*drkp1*dt;

    % Jacobian
    dR = I - (dfxkp1 + dfzkp1*dzkp1)*dt;
else
    % Right-hand side function
    fkp1 = f(tkp1, xkp1, zkp1, p);
end

% Evaluate residual equations
R = xkp1 - xk - fkp1*dt;

% Rescale residual function
R = R/dt;
dR = dR/dt;