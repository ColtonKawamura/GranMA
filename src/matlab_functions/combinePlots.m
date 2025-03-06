function combinePlots(f1, f2)
    % Parameters:
    % f1 - Filename of the first .fig file in "string" format
    % f2 - Filename of the second .fig file

    % Define the two figure files that are on the MATLAB path
    % These are the only "inputs" to this script

    % Open the figures
    fig1 = openfig(f1);
    fig2 = openfig(f2);

    % get axes handles - this assumes there is only 1 axes per figure!
    fig1ax = gca(fig1);
    fig2ax = gca(fig2);
    leg1 = findobj(fig1,'Type','legend');
    % leg2 = findobj(fig2,'Type','legend');

    % Get axis children
    fig1axChildren = get(fig1ax,'Children');
    fig2axChildren = get(fig2ax,'Children');

    % Create new fig, copy items from fig 1
    % This will maintain set properties such as color
    figFinal = figure();
    ax = axes(figFinal);
    h1 = copyobj(fig1axChildren, ax);

    % Copy items from fig 2
    h2 = copyobj(fig2axChildren, ax);
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    grid on

    % change the maker shapes

    set(h1, 'Marker', 'o', 'MarkerFaceColor', 'none');
    set(h2, 'Marker', '^', 'MarkerFaceColor', 'none');
    % % Add legend to same location as the legend in fig2 
    % % but only include objects with a defined DisplayName
    h = [h2;h1];
    hasDisplayName = ~cellfun('isempty',get(h,'DisplayName'));
    leg = legend(ax, h(hasDisplayName),'Location', leg1.Location, 'Interpreter', 'Latex')
    % title(leg, " $\hat{P} $")

    % Copy axis labels
    xlabel(ax, fig2ax.XLabel.String, 'Interpreter', 'Latex', 'FontSize', 20)
    ylabel(ax, fig2ax.YLabel.String, 'Interpreter', 'Latex', 'FontSize', 20)
    set(get(gca, 'ylabel'), 'rotation', 0);
    box on;

    % Legend from a single simulation (pulls from figure_2 only)
    % fig_handle = get(fig2, 'Children');
    % hasDisplayName = ~cellfun('isempty', get(fig_handle, 'DisplayName'));
    % legend_entries = fig_handle(hasDisplayName);
    % leg = legend(fig, legend_entries, 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeastoutside', 'Orientation', 'vertical');
    % title(leg, "$ \hat{P}$ ")
    % Close original figures
    close(fig1);
    close(fig2);

end
