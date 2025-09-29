function lambda = steady_state_eigenvalues_subset_ddde(x, f, h, p, c, a, M, varargin)
% Evaluate the delayed quantity
[r, dhx] = h(x, p, varargin{:});

% Evaluate the delayed state
z = r;

% Variable dimensions
nx = numel(x);
nz = numel(z);

% Evaluate the residual function
[~, dfx, dfz] = f([], x, z, p, varargin{:});

% Dimension of matrix
dimA = nx + sum(M) + nz;

% Function handle
Afun = @(x) A(x, dfx, dfz, dhx, c, a, M);

            % Abar = [dfx, dfz*p.C; p.B*dhx, p.A];
            % x = randn(size(Abar, 1), 1);
            % if(norm(Abar*x - Afun(x))/norm(Abar*x) > 1e2*eps/2), error('!'); end

% Compute a subset of the eigenvalues
lambda = eigs(Afun, dimA, varargin{:});
end

function Ainp = A(inp, dfx, dfz, dhx, c, a, M)
% Variable dimensions
nx = size(dfx, 2);
nz = size(dfz, 2);

% Auxiliary variables used to identify indices
aux = [nx; M+1];
idx = [[1; cumsum(aux(1:end-1))+1], cumsum(aux)];

% Local indices
idx_start = idx(1, 1);
idx_end   = idx(1, 2);

% First row block of input
inp_0 = inp(idx_start:idx_end);

% Precomputation
dhx_x_0 = dhx*inp_0;

% Allocate memory
C_inp_i = zeros(        nz, 1);
Ainp    = zeros(numel(inp), 1);

for i = 1:numel(M)
    % Local indices
    idx_start = idx(i+1, 1);
    idx_end   = idx(i+1, 2);

    % Local vector
    inp_i = inp(idx_start:idx_end);

    % Precomputation used below
    C_inp_i(i, 1) = c(i, 1:M(i)+1)*inp_i;

    % Insert expression
    Ainp(idx_start:idx_end, 1) = a(i)*([dhx_x_0(i); inp_i(1:end-1)] - inp_i);
end

% Local indices
idx_start = idx(1, 1);
idx_end   = idx(1, 2);

% First row block
Ainp(idx_start:idx_end, 1) = dfx*inp_0 + dfz*C_inp_i;
end