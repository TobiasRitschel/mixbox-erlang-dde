function SavePlot(name,varargin)
% Save plot in a given format
%
% SYNOPSIS:
%   SavePlot(name,varargin)
%
% DESCRIPTION:
% Saves figure in eps (default), tikz or png format. eps figures are saved
% in a subdirectory called eps, tikz are saved in a subdirectory called
% tikz and png are saved in a subdirectory called png. If the relevant
% subdirectory does not exist SavePlot will generate an error.
%
% REQUIRED PARAMETERS:
%   name        - name of saved file
%   varargin    - 1: boolean indicating whether to save or not
%                 2: format: eps (default), tikz or png
%
% RETURNS:

% Default settings
save = false;
format = 'eps';

% Check if user wants to save the plot
if(nargin > 1)
    save = varargin{1};
end

% Check if the user has specified the format
if(nargin > 2)
    format = varargin{2};
end

% Save if necessary
if(save)
    switch lower(format)
        case 'tikz'
            matlab2tikz(['tikz/' name '.tikz'], ...
                'height','\figureheight','width','width','\figurewidth', ...
                'ShowInfo',false, 'ShowWarnings',false);
        case 'eps'
            saveas(gcf,['eps/' name],'epsc');
        case 'png'
            print(['png/' name],'-dpng','-r600');
        case 'pdf'
            % Source: https://stackoverflow.com/questions/3801730/get-rid-of-the-white-space-around-matlab-figures-pdf-output
            h = gcf;
            set(h, 'PaperUnits','centimeters');
            set(h, 'Units','centimeters');
            pos=get(h,'Position');
            set(h, 'PaperSize', [pos(3) pos(4)]);
            set(h, 'PaperPositionMode', 'manual');
            set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
            print('-dpdf', ['pdf/' name]);
%             saveas(gcf,['pdf/' name, '.pdf']);
    end
end