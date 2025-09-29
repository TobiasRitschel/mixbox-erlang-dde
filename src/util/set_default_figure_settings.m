function set_default_figure_settings(fs, lw, ms, fig_pos)
if(~isunix)
    % Set default font size
    set(groot, 'DefaultAxesFontSize',   fs);

    % Set default line widths
    set(groot, 'DefaultLineLineWidth',  lw);
    set(groot, 'DefaultStairLineWidth', lw);
    set(groot, 'DefaultStemLineWidth',  lw);

    % Set default marker size
    set(groot, 'DefaultLineMarkerSize', ms);

    % Set default interpreter
    set(groot, 'DefaultTextInterpreter',            'Latex');
    set(groot, 'DefaultAxesTickLabelInterpreter',   'Latex');
    set(groot, 'DefaultLegendInterpreter',          'Latex');

    % Set default figure size
    set(groot, 'DefaultFigureUnits', 'Inches', 'DefaultFigurePosition', fig_pos);

    % Set default renderer (prevent pixelated figures from eps)
    set(groot, 'DefaultFigureRenderer', 'Painters');
end