function beta = evaluateKernelIntegral(t, c, a)

% Highest polynomial order
Mp1 = numel(c);

% Initialize
beta = 0;

for idx = 1:Mp1
    % Local index
    m = idx-1;

    % Create auxiliary weights
    c_inner = ones(m+1, 1);

    % Evaluate auxiliary kernel
    alpha_inner = evaluateKernel(t, c_inner, a);
    
    % Compute m'th integral
    betam = 1 - alpha_inner/a;

    % Add to integral
    beta = beta + c(idx)*betam;
end