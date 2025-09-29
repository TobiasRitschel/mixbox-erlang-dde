function [out1, out2] = steady_state_eigenvalues_ddde(x, f, h, p, varargin)
% Evaluate the delayed quantity
[r, dhx] = h(x, p, varargin{:});

% Evaluate the delayed state
z = r;

% Evaluate the residual function
[~, dfx, dfz] = f([], x, z, p, varargin{:});

% Jacobian
Abar = [dfx, dfz*p.C; p.B*dhx, p.A];

if(nargout > 1)
    % Compute all eigenvalues and eigenvectors
    [V, D] = eig(Abar);

    % Extract eigenvalues
    lambda = diag(D);

    % Assign outputs
    out1 = V;
    out2 = lambda;
else
    % Compute all eigenvalues and eigenvectors
    lambda = eig(Abar);

    % Assign outputs
    out1 = lambda;
end