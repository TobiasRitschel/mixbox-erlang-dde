function [alpha, dalpha, d2alpha, alpham_out] = evaluateKernel(t, c, a)

% Compute gradient and Hessian?
ComputeGradient = (nargout > 1);
ComputeHessian  = (nargout > 2);

%% Mixed Erlang kernel
% Number of kernels
nz = numel(a);

% Number of time points
N = numel(t);

% Highest polynomial order
Mp1 = size(c, 2);

% Initialize
alpha = zeros(nz, N);

if(ComputeGradient)
    % Allocate memory
    dalpha = zeros(nz, numel(t), Mp1+1);

    if(ComputeHessian)
        % Allocate memory
        d2alpha = zeros(nz, numel(t), Mp1+1, Mp1+1);
    end
end

% Initialize
bm = 1;

for idx = 1:Mp1
    % Local index
    m = idx-1;

    % Normalization factor
    bm = a.*bm/max(1, m);

    % Erlang probability density function (m'th order)
    alpham = bm.*t.^m.*exp(-a*t);

    % Kernel
    alpha = alpha + c(:, idx).*alpham;

    if(nargout > 3)
        % Store subkernel
        alpham_out(idx, :) = alpham; %#ok
    end

    if(ComputeGradient)
        % Derivative of probability density function
        dalpham_a = ((m+1)./a - t).*alpham;

        % Derivative of kernel wrt. coefficient
        dalpha(:, :, idx) = alpham;

        % Derivative of kernel wrt. rate parameter
        dalpha(:, :, Mp1+1) = dalpha(:, :, Mp1+1) + c(:, idx).*dalpham_a;

        if(ComputeHessian)
            % Derivative of probability density function
            d2alpham_a = -(m+1)./a.^2.*alpham + ((m+1)./a - t).*dalpham_a;

            % Mixed derivative
            d2alpha(:, :, idx, Mp1+1) = dalpham_a;
            d2alpha(:, :, Mp1+1, idx) = dalpham_a;

            % Second order derivative of kernel wrt. rate parameter
            d2alpha(:, :, Mp1+1, Mp1+1) = squeeze(d2alpha(:, :, Mp1+1, Mp1+1)) + c(:, idx).*d2alpham_a;
        end
    end
end