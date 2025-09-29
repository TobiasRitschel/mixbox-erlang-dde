function hess = least_squares_kernel_objective_hessian(theta, ~, t, alpha, opts)

% Number of parameters
ntheta = numel(theta);

% Extract parameters
c = reshape(theta(1:end-1), 1, []);
a = reshape(theta(  end  ), 1, []);

% Scaling factors
if(~isempty(opts) && isfield(opts, 'rho')), rho = opts.rho; else, rho = 1; end

% Indices of points used in quadrature
idx = 1:numel(t)-1;

% Choose the left points in the integral discretization
tl = t(idx);

% Compute the time intervals
dt = diff(t);

% Evaluate approximate kernel
[alphahat, dalphahat, d2alphahat] = evaluateKernel(tl, c, a);

% Compute deviation
e = alpha(idx) - alphahat;

% Derivative of deviation
de = -dalphahat;

% Hessian of deviation
d2e = -d2alphahat;

% Allocate memory
hess = zeros(ntheta, ntheta);

for i = 1:ntheta
    for j = 1:ntheta
        % Hessian element
        hess(i, j) = rho*((de(:, :, j).*dt)*de(:, :, i)' + (e.*dt)*d2e(:, :, i, j)');
    end
end