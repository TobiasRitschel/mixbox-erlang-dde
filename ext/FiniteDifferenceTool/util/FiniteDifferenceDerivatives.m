function [dfFD, d2fFD] = FiniteDifferenceDerivatives(fun, x, epsilon, varargin)
% FiniteDifferenceDerivatives  Currently undocumented

% Evaluate function
f = feval(fun, x, varargin{:});

% Dimensions
nx = length(x);
nf = numel(f);

%% First order derivatives
e = zeros(nx, 1);
dfFD = zeros(nf, nx);
for j = 1:nx
    % Pertubation
    e(j) = epsilon;
    
    % Perturbed function evaluation
    fFDp = feval(fun, x+e, varargin{:});
    fFDm = feval(fun, x-e, varargin{:});
    
    % Approximation
    for i = 1:nf
        dfFD(i,j) = (fFDp(i) - fFDm(i)) / (2*epsilon);
    end
    
    % Remove pertubation
    e(j) = 0;
end

%% Second order derivatives
if(nargout > 1)
    d2fFD = zeros(nx, nx, nf);
    for j = 1:nx
        for k = 1:nx
            % Pertubation
            e = zeros(nx, 1);
            e(j) =        epsilon;
            e(k) = e(k) + epsilon;
            
            % Perturbed function evaluation
            fFDpp = feval(fun, x+e, varargin{:});
            fFDmm = feval(fun, x-e, varargin{:});
            
            % Pertubation
            e = zeros(nx, 1);
            e(j) =        epsilon;
            e(k) = e(k) - epsilon;
            
            % Perturbed function evaluation
            fFDpm = feval(fun, x+e, varargin{:});
            
            % Pertubation
            e = zeros(nx, 1);
            e(j) =      - epsilon;
            e(k) = e(k) + epsilon;
             
            % Perturbed function evaluation
            fFDmp = feval(fun, x+e, varargin{:});
            
            % Approximation
            for i = 1:nf
                d2fFD(j, k, i) = (fFDpp(i) - fFDpm(i) - fFDmp(i) + fFDmm(i)) / (4*epsilon*epsilon);
            end
        end
    end
end