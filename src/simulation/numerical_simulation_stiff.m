function [T, X, Z, R] = numerical_simulation_stiff(f, tspan, x0, p, alpha, h, Nh, opts)

% Initial delayed quantity (overwritten below)
r0 = h(x0, p);

% Number of states and delayed states
nx = size(x0, 1);
nz = size(r0, 1);

% Simulation horizon
N = numel(tspan)-1;

% Time step size
dt = diff(tspan(1:2)); % Note: tspan must be a vector of equidistant times

% Evaluate kernel
alphatkmsj = alpha((0:Nh-1)*dt);

% Times
T = (-Nh+1:N)*dt;

% Allocate memory
X = zeros(nx, N+Nh);
Z = zeros(nz, N+Nh);
R = zeros(nz, N+Nh);

if(size(x0, 2) > 1)
    % Insert initial state
    X(:, 1:Nh) = x0;

    for k = 1:Nh
        % Initial delayed quantity
        r0 = h(x0(:, k), p);

        % Insert initial delayed quantity
        R(:, k) = r0;
    end

    % Initialize
    xk = x0(:, end);
else
    for k = 1:Nh
        % Insert initial condition
        X(:, k) = x0;
        R(:, k) = r0;
    end

    % Initialize
    xk = x0;
end

for k = Nh + (0:N-1)
    % Time
    tkp1 = T(k+1);

    % Historical states (excluding rk+1)
    rkmjp1 = R(:, k - (1:Nh-1) + 1);

    % Compute states by solving the residual equations
    xkp1 = fsolve(@residual_equations, xk, opts, ...
        tkp1, xk, rkmjp1, alphatkmsj, dt, p, f, h);

    % Delayed quantities
    rkp1 = h(xkp1, p);

    % Compute delayed states
    zkp1 = sum(alphatkmsj.*[rkp1, rkmjp1], 2)*dt;

    % Store solution
    X(:, k+1) = xkp1;
    Z(:, k+1) = zkp1;
    R(:, k+1) = rkp1;

    % Update initial values
    xk = xkp1;
end
end