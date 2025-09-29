function [T, X, Z, R] = numerical_simulation_nonstiff(f, tspan, x0, p, alpha, h, Nh)

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
alphatkp1msj = alpha((1:Nh)*dt);

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

    % Historical states
    r0mj = R(:, Nh - (1:Nh-1));

    % Initialize
    xk = x0(:, end);
    zk = sum(alphatkp1msj(1:end-1).*r0mj)*dt;
else
    for k = 1:Nh
        % Insert initial condition
        X(:, k) = x0;
        R(:, k) = r0;
    end

    % Initialize
    xk = x0;
    zk = r0;
end

for k = Nh + (0:N-1)
    % Time
    tk = T(k);

    % Historical states
    rkmjp1 = R(:, k - (1:Nh) + 1);

    % Compute delayed states
    zkp1 = sum(alphatkp1msj.*rkmjp1, 2)*dt;

    % Compute states
    xkp1 = xk + f(tk, xk, zk, p)*dt;

    % Store solution
    X(:, k+1) = xkp1;
    Z(:, k+1) = zkp1;
    R(:, k+1) = h(xkp1, p);

    % Update initial values
    xk = xkp1;
    zk = zkp1;
end