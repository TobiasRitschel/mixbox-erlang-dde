function stable = steady_state_stability_ddde(x, f, h, p, varargin)

% Evaluate the delayed quantity
[r, dhx] = h(x, p, varargin{:});

% Evaluate the delayed state
z = r;

% Evaluate the residual function
[~, dfx, dfz] = f([], x, z, p, varargin{:});

% Jacobian
Abar = [dfx, dfz*p.C; p.B*dhx, p.A];

% Compute the eigenvalues
lambda_max_real_part = eigs(Abar, 1, 'largestreal');

% Assess stability
stable = lambda_max_real_part < 0;
end