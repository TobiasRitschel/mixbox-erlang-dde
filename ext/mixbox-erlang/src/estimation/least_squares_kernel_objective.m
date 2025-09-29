function [phi, grad] = least_squares_kernel_objective(theta, t, alpha, opts)
% Extract parameters
c = reshape(theta(1:end-1), 1, []);
a = reshape(theta(  end  ), 1, []);

% Scaling factors
if(~isempty(opts) && isfield(opts, 'rho')), rho = opts.rho; else, rho = 1; end

% Number of measurements
N = numel(t);

% Indices of points used in quadrature
idx = 1:N-1;

% Choose the left points in the integral discretization
tl = t(idx);

% Compute the time intervals
dt = diff(t);

if(nargout > 1)
    % Evaluate approximate kernel
    [alphahat, dalphahat] = evaluateKernel(tl, c, a);
else
    % Evaluate approximate kernel
    alphahat = evaluateKernel(tl, c, a);
end

% Compute deviation
e = alpha(idx) - alphahat;

% Precomputation
edt = e.*dt;

% Evaluate objective function
phi = 0.5*rho*(edt*e');

if(nargout > 1)
    % Derivative of deviation
    de = -dalphahat;

    % Derivative of objective function
    grad = rho*edt*reshape(permute(de, [2, 3, 1]), N-1, []);
end