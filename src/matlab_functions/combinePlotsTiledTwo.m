function combinePlotsTiledTwo(f1, f2, x_limits, y_limits)
    % combinePlotsTiledTwo Combines two .fig files in a tiled layout with a single legend.
    % 
    % Parameters:
    % f1 - Filename of the first .fig file in "string" format
    % f2 - Filename of the second .fig file
    % x_limits = [lower, upper] (Optional)
    % y_limits = [lower, upper] (Optional)

    if nargin < 3 || isempty(x_limits)
        x_limits = [];  % Default: No limits set
    end
    if nargin < 4 || isempty(y_limits)
        y_limits = [];  % Default: No limits set
    end
    
    % Load figures
    fig1 = openfig(f1, 'invisible');
    fig2 = openfig(f2, 'invisible');

    % Get axes and children for each figure
    fig1ax = gca(fig1);
    fig2ax = gca(fig2);
    fig1axChildren = get(fig1ax, 'Children');
    fig2axChildren = get(fig2ax, 'Children');

    % Create the main figure with tiled layout
    figure_main = figure;
    tiled_main = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'none');

    % First tile - Plot 1
    figure_1 = nexttile;
    copyobj(fig1axChildren, figure_1);
    grid(figure_1, 'on');
    box(figure_1, 'on');  
    ylabel(figure_1, fig1ax.YLabel.String, 'Interpreter', 'latex', 'FontSize', 15);
    set(figure_1, 'XTickLabel', ''); % Hide x-tick labels but keep tick marks visible
    if ~isempty(x_limits)
        xlim(figure_1, x_limits);
    end

    % Remove y-tick at 0
    y_ticks = yticks();
    y_ticks(y_ticks == 0) = [];
    yticks(y_ticks);

    % Second tile - Plot 2
    figure_2 = nexttile;
    copyobj(fig2axChildren, figure_2);
    grid(figure_2, 'on');
    box(figure_2, 'on');  
    set(figure_2, 'YScale', 'log'); % Log scale for y-axis
    ylabel(figure_2, fig2ax.YLabel.String, 'Interpreter', 'latex', 'FontSize', 15);
    xlabel(figure_2, fig2ax.XLabel.String, 'Interpreter', 'latex', 'FontSize', 15);
    
    % Match x-limits with the first plot
    xLimits = xlim(figure_1);
    xlim(figure_2, xLimits);
    if ~isempty(y_limits)
        ylim(figure_2, y_limits);
    end

    % Remove y-tick at 0
    y_ticks = yticks();
    y_ticks(y_ticks == 0) = [];
    yticks(y_ticks);
    
    % Extract legend from figure_2
    fig_handle = get(figure_2, 'Children');
    hasDisplayName = ~cellfun('isempty', get(fig_handle, 'DisplayName'));
    legend_entries = fig_handle(hasDisplayName);
    legend(figure_2, legend_entries, 'Interpreter', 'latex', 'FontSize', 12, ...
        'Location', 'northeastoutside', 'Orientation', 'vertical');

    % Close original figures
    close(fig1);
    close(fig2);
end
