function th = identify_domain_bisection_recursion(beta, epsilon, tl, tu, tol, varargin)
% Compute midpoint
tm = 0.5*(tu + tl);

% Evaluate cumulative distribution function at the midpoint
Fm = beta(tm, varargin{:});

% Check for convergence
Converged = abs(1 - Fm - epsilon) < tol;

% Return if the algorithm has converged
if(Converged), th = tm; return; end

if(1 - Fm < epsilon)
    % Repeat for the lower interval
    th = identify_domain_bisection_recursion(beta, epsilon, tl, tm, tol, varargin{:});
else
    % Repeat for the upper interval
    th = identify_domain_bisection_recursion(beta, epsilon, tm, tu, tol, varargin{:});
end
end