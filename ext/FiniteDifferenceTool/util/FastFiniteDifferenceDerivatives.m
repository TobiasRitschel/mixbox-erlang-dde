function [dfFD, d2fFD] = FastFiniteDifferenceDerivatives(fun, x, epsilon, varargin)
% FiniteDifferenceDerivatives  Currently undocumented

% Evaluate function
f = feval(fun, x, varargin{:});

% Dimensions
nx = numel(x);
nf = numel(f);

%% First order derivatives
dfFD = zeros(nf, nx);
for j = 1:nx
    % Pertubation
    x(j) = x(j) + epsilon;
    
    % Perturbed function evaluation
    fp = feval(fun, x, varargin{:});
    
    % Approximation
    dfFD(:, j) = (fp(:) - f(:)) / epsilon;
    
    % Remove pertubation
    x(j) = x(j) - epsilon;
end

%% Second order derivatives
if(nargout > 1)
    epssq = epsilon*epsilon;
    d2fFD = zeros(nx, nx, nf);
    for j = 1:nx
        % Pertubation
        x(j) = x(j) + 2*epsilon;
            
        % Perturbed function evaluation
        fpp = feval(fun, x, varargin{:});
            
        % Pertubation
        x(j) = x(j) - epsilon;

        % Perturbed function evaluation
        fpz = feval(fun, x, varargin{:});
        
        % Approximation
        d2fFD(j, j, :) = (fpp(:) - 2*fpz(:) + f(:)) / epssq;
        
        % Reset pertubation
        x(j) = x(j) - epsilon;
        
        for k = 1:j-1
            % Pertubation
            x(j) = x(j) + epsilon;
            x(k) = x(k) + epsilon;
            
            % Perturbed function evaluation
            fpp = feval(fun, x, varargin{:});
            
            % Reset pertubation
            x(k) = x(k) - epsilon;
             
            % Perturbed function evaluation
            fpz = feval(fun, x, varargin{:});
            
            % Pertubation
            x(k) = x(k) + epsilon;
            x(j) = x(j) - epsilon;
             
            % Perturbed function evaluation
            fzp = feval(fun, x, varargin{:});
            
            % Approximation
            d2fFD(j, k, :) = (fpp(:) - fpz(:) - fzp(:) + f(:)) / epssq;
            d2fFD(k, j, :) = d2fFD(j, k, :);
            
            % Reset pertubation
            x(k) = x(k) - epsilon;
        end
    end
end