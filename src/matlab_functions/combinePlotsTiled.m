function combinePlotsTiled(f1, f2)
    % combinePlotsTiled Combines two .fig files in a tiled layout with a single legend.
    % 
    % Parameters:
    % f1 - Filename of the first .fig file in "string" format
    % f2 - Filename of the second .fig file

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

    % First tile - general plot 1
    figure_1 = nexttile;
    copyobj(fig1axChildren, figure_1);
    set(figure_1, 'YScale', 'log', 'XScale', 'log');
    grid(figure_1, 'on');
    box(figure_1, 'on'); % Turn on the box for the first tile
    ylabel(figure_1, fig1ax.YLabel.String, 'Interpreter', 'latex', 'FontSize', 15); % Original y-label for fig1
    set(get(figure_1, 'ylabel'), 'rotation', 0);
    set(figure_1, 'XTickLabel', ''); % Hide x-tick labels but keep tick marks visible
    xlim(figure_1, [.013, 3])

    % Second tile - general plot 2
    figure_2 = nexttile;
    copyobj(fig2axChildren, figure_2);
    grid(figure_2, 'on');
    box(figure_2, 'on'); % Turn on the box for the second tile
    set(figure_2, 'YScale', 'log', 'XScale', 'log'); % Both y and x log
    set(figure_2, 'XScale', 'log'); % just x asix
    ylabel(figure_2, fig2ax.YLabel.String, 'Interpreter', 'latex', 'FontSize', 15); % Original y-label for fig2
    xlabel(figure_2, fig2ax.XLabel.String, 'Interpreter', 'latex', 'FontSize', 15); % Original x-label for fig2
    % set(get(figure_2, 'ylabel'), 'rotation', 0);

    % Match x-limits
    ylim(figure_1, [.03, .7])
    xLimits = xlim(figure_1);  % Get x-limits from the first plot
    xlim(figure_2, xLimits);   % Apply the same x-limits to the second plot

    % Legend from a single simulation (pulls from figure_2 only)
    fig_handle = get(figure_2, 'Children');
    hasDisplayName = ~cellfun('isempty', get(fig_handle, 'DisplayName'));
    legend_entries = fig_handle(hasDisplayName);
    leg = legend(figure_2, legend_entries, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeastoutside', 'Orientation', 'vertical');
    title(leg, "$ \hat{P}$ ")
    % Close original figures
    close(fig1);
    close(fig2);
end
