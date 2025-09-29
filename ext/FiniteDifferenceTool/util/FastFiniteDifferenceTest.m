function [edf, ed2f] = FastFiniteDifferenceTest(fun, x, epsilon, PlottingOn, ErrorFunction, varargin)
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
nx = numel(x);
nf = numel(f);

% Finite differences
NoTests = length(epsilon);

% Finite difference tests
etotdf  = zeros(NoTests, 1);
edf     = zeros(NoTests, nx*nf);
if(nout > 2)
    etotd2f = zeros(NoTests, 1);
    ed2f    = zeros(NoTests, nx*nx*nf);
end
tic
for i = 1:NoTests
    % Approximate derivatives
    if(nout > 2)
        [dfFD, d2fFD] = FastFiniteDifferenceDerivatives(fun, x(:), epsilon(i), varargin{:});
    else
        dfFD = FastFiniteDifferenceDerivatives(fun, x(:), epsilon(i), varargin{:});
    end
    
    % Evaluate first order approximation
    edf(i, :) = ErrorFunction(df(:), dfFD(:));
    
    % Detect zero derivatives
    idx = abs(dfFD(:)) > sqrt(eps);
    
    % Sum error of nonzero derivatives
    etotdf(i) = norm(edf(i,idx));
    
    if(nout > 2)
        % Evaluate first order approximation
        ed2f(i,:) = ErrorFunction(d2f(:), d2fFD(:));
        
        idx = abs(d2fFD(:)) > sqrt(eps);
        etotd2f(i) = norm(ed2f(i,idx));
    end
    
    if(i > 1), fprintf(repmat('\b', 1, 14)); end
    fprintf('Progress: %3.0f%%', i/NoTests*1e2);
end
fprintf('\n');

%% Visualize results
if(PlottingOn)
    if(nx*nf <= 5e2)
        figure;
        clf;
        loglog(epsilon, edf, '-o');
        xlabel('\epsilon')
        ylabel('(x - xFD)./xFD');
        title('Individual 1st order errors');
        legend(num2str((1:nx*nf)'), 'location', 'northwest');
        
        if(nout > 2)
            figure;
            clf;
            loglog(epsilon, ed2f, '-x');
            xlabel('\epsilon')
            ylabel('(x - xFD)./xFD');
            title('Individual 2nd order errors');
            legend(num2str((1:nx*nx*nf)'), 'location', 'northwest');
        end
    end
    
    figure;
    clf;
    loglog(epsilon, etotdf, '-bx', 'DisplayName', '1st der.'); hold on
    if(nout > 2), loglog(epsilon, etotd2f, '-ro', 'DisplayName', '2nd der.'); hold off; end
    xlabel('\epsilon')
    ylabel('||(x - xFD)./xFD||');
    legend('location','northwest');
    title('Total errors');
end