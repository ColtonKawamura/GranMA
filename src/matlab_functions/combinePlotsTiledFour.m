function combinePlotsTiledFour(f1, f2, f3, f4, x_limits, y_limits)
    % combinePlotsTiled Combines two .fig files in a tiled layout with a single legend.
    % 
    % Parameters:
    % f1 - Filename of the first .fig file in "string" format
    % f2 - Filename of the second .fig file
    % etc.
    % xlim = [lower, uppper]
    % ylim = [lower, uppper]

    if nargin < 5 || isempty(x_limits)
        x_limits = [];  % Default: No limits set
    end
    if nargin < 6 || isempty(y_limits)
        y_limits = [];  % Default: No limits set
    end
    
    % Load figures
    fig1 = openfig(f1, 'invisible');
    fig2 = openfig(f2, 'invisible');
    fig3 = openfig(f3, 'invisible');
    fig4 = openfig(f4, 'invisible');

    % Get axes and children for each figure
    fig1ax = gca(fig1);
    fig2ax = gca(fig2);
    fig3ax = gca(fig3);
    fig4ax = gca(fig4);
    fig1axChildren = get(fig1ax, 'Children');
    fig2axChildren = get(fig2ax, 'Children');
    fig3axChildren = get(fig3ax, 'Children');
    fig4axChildren = get(fig4ax, 'Children');

    % Get x and y limits from both figures
    x1Limits = xlim(fig1ax);
    y1Limits = ylim(fig1ax);
    x2Limits = xlim(fig2ax);
    y2Limits = ylim(fig2ax);
    x3Limits = xlim(fig3ax);
    y3Limits = ylim(fig3ax);
    x4Limits = xlim(fig4ax);
    y4Limits = ylim(fig4ax);

    % Determine the global min/max for x and y
    xMin = min(x1Limits(1), x2Limits(1));
    xMax = max(x1Limits(2), x2Limits(2));
    yMin = min(y1Limits(1), y2Limits(1));
    yMax = max(y1Limits(2), y2Limits(2))

    % Create the main figure with tiled layout
    figure_main = figure;
    tiled_main = tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'none');

    % First tile - general plot 1
    figure_1 = nexttile;
    copyobj(fig1axChildren, figure_1);
    % set(figure_1, 'YScale', 'log', 'XScale', 'log');
    grid(figure_1, 'on');
    box(figure_1, 'on'); % Turn on the box for the first tile
    ylabel(figure_1, fig1ax.YLabel.String, 'Interpreter', 'latex', 'FontSize', 15); % Original y-label for fig1
    set(get(figure_1, 'ylabel'), 'rotation', 0);
    set(figure_1, 'XTickLabel', ''); % Hide x-tick labels but keep tick marks visible
    if ~isempty(x_limits)
        xlim(figure_1, x_limits) %
    end
    % if ~isempty(y_limits)
    %     xlim(figure_1, y_limits) %
    % end
    % xlim(figure_1, [.001, 1.1]) % for gamma = .01
    % ylim(figure_1, [.01, 1.5]);
        % Apply global x and y limits to the first plot
    % xlim(figure_1, [xMin, xMax]);
    % ylim(figure_1, [yMin, yMax]);


    % Second tile - general plot 2
    figure_2 = nexttile;
    copyobj(fig2axChildren, figure_2);
    grid(figure_2, 'on');
    box(figure_2, 'on'); % Turn on the box for the second tile
    % set(figure_2, 'YScale', 'log', 'XScale', 'log'); % Both y and x log
    set(figure_2, 'YScale', 'log'); % just x asix
    ylabel(figure_2, fig2ax.YLabel.String, 'Interpreter', 'latex', 'FontSize', 15); % Original y-label for fig2
    set(figure_2, 'XTickLabel', ''); % Hide x-tick labels but keep tick marks visiblej

    % Match x-limits
    % ylim(figure_1, [.03, .7])
    xLimits = xlim(figure_1);  % Get x-limits from the first plot
    xlim(figure_2, xLimits);   % Apply the same x-limits to the second plot
    if ~isempty(y_limits)
        ylim(figure_2, y_limits)
    end
   
    
    % Third tile - general plot 3
    figure_3 = nexttile;
    copyobj(fig3axChildren, figure_3);
    grid(figure_3, 'on');
    box(figure_3, 'on'); % Turn on the box for the second tile
    % set(figure_2, 'YScale', 'log', 'XScale', 'log'); % Both y and x log
    % set(figure_3, 'YScale', 'log'); % just x asix
    ylabel(figure_3, fig3ax.YLabel.String, 'Interpreter', 'latex', 'FontSize', 15); % Original y-label for fig2
    set(figure_3, 'XTickLabel', ''); % Hide x-tick labels but keep tick marks visible

    % Match x-limits
    % ylim(figure_1, [.03, .7])
    xLimits = xlim(figure_1);  % Get x-limits from the first plot
    xlim(figure_3, xLimits);   % Apply the same x-limits to the second plot
    % if ~isempty(y_limits)
    %     ylim(figure_3, y_limits)
    % end

    % Fourth tile - general plot 4
    figure_4 = nexttile;
    copyobj(fig4axChildren, figure_4);
    grid(figure_4, 'on');
    box(figure_4, 'on'); % Turn on the box for the second tile
    % set(figure_2, 'YScale', 'log', 'XScale', 'log'); % Both y and x log
    set(figure_4, 'YScale', 'log'); % just x asix
    ylabel(figure_4, fig4ax.YLabel.String, 'Interpreter', 'latex', 'FontSize', 15); % Original y-label for fig2
    xlabel(figure_4, fig4ax.XLabel.String, 'Interpreter', 'latex', 'FontSize', 15); % Original x-label for fig2
    % set(get(figure_2, 'ylabel'), 'rotation', 0);


    % Match x-limits
    % ylim(figure_1, [.03, .7])
    xLimits = xlim(figure_1);  % Get x-limits from the first plot
    xlim(figure_4, xLimits);   % Apply the same x-limits to the second plot
    if ~isempty(y_limits)
        ylim(figure_4, y_limits)
    end

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
