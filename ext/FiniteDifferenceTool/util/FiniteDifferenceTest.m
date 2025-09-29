function [edf, ed2f] = FiniteDifferenceTest(fun, x, epsilon, PlottingOn, ErrorFunction, varargin)
% FiniteDifferenceTest  Currently undocumented

% Check if error function is defined
if(isempty(ErrorFunction))
    ErrorFunction = @(df, dfFD) abs((df(:) - dfFD(:))./dfFD(:));
end

% Number of output arguments
nout = nargout(fun);

% Evalute function
if(nout > 2)
    [f, df, d2f] = fun(x, varargin{:});
else
    [f, df] = fun(x, varargin{:});
end

% Dimensions
nx = length(x);
nf = numel (f);

% Finite differences
NoTests = length(epsilon);

% Finite difference tests
etotdf = zeros(NoTests, 1);
etotd2f = zeros(NoTests, 1);
edf = zeros(NoTests, nx*nf);
ed2f = zeros(NoTests, nx*nx*nf);
for i = 1:NoTests
    % Approximate derivatives
    [dfFD, d2fFD] = FiniteDifferenceDerivatives(fun, x, epsilon(i), varargin{:});

    % Evaluate first order approximation
    edf(i,:) = ErrorFunction(df(:), dfFD(:));
    
    idx = abs(df(:)) > sqrt(eps);
    etotdf(i) = norm(edf(i,idx));
    
    if(nout > 2)
        % Evaluate second order approximation
        ed2f(i,:) = ErrorFunction(d2f(:), d2fFD(:));
        
        idx = abs(d2f(:)) > sqrt(eps);
        etotd2f(i) = norm(ed2f(i,idx));
    end
    
    if(i > 1), fprintf(repmat('\b', 1, 14)); end
    fprintf('Progress: %3.0f%%', i/NoTests*1e2);
end
fprintf('\n');

%% Visualize results
if(PlottingOn)
    figure;
    clf;
    loglog(epsilon, edf, '-o');
    xlabel('\epsilon')
    ylabel('Error');
    title('Individual 1st order errors');
    legend(num2str((1:nx*nf)'), 'location', 'northwest');
    
    if(nout > 2)
        figure;
        clf;
        loglog(epsilon, ed2f, '-x');
        xlabel('\epsilon')
        ylabel('Error');
        title('Individual 2nd order errors');
        legend(num2str((1:nx*nx*nf)'), 'location', 'northwest');
    end
    
    figure;
    clf;
    loglog(epsilon, etotdf, '-bx');
    if(nout > 2), hold on; loglog(epsilon, etotd2f, '-ro'); hold off; end
    xlabel('\epsilon')
    ylabel('Error');
    if(nout > 2), legend('1st der.', '2nd der.','location','northwest'); end
    title('Total errors');
end