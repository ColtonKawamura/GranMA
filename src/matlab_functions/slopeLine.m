function slopeLine(type, slope, x_bounds, y_center, options)
    % SLOPELINE Plot a line with a given slope and center point, with text label placement.
    %   SLOPELINE(TYPE, SLOPE, X_BOUNDS, Y_CENTER) plots a line with default text placement.
    %   SLOPELINE(..., 'TextLocation', [X, Y]) specifies the exact coordinates for the text label.

    arguments
        type {mustBeMember(type, {'linear', 'semilog', 'loglog'})}
        slope (1,1) double
        x_bounds (1,:) double % Ensure x_bounds is a row vector
        y_center (1,1) double
        options.TextLocation (1,2) double = [] % Optional: [x,y] coordinates for text
    end

    % ... (Input validation and line calculation code remains the same) ...
    % Ensure x_bounds has at least two elements for min/max/mean
    if numel(x_bounds) < 2
        error('x_bounds must contain at least two elements.');
    end
    % Ensure x_bounds is sorted for correct y(1) and y(end) mapping
    x = sort(x_bounds);

    hold on;
    % Calculate y values based on the type of plot
    switch type
        case 'linear'
            y = slope * (x - mean(x)) + y_center;
        case 'semilog'
            % Ensure x values are positive for log10
            if any(x <= 0)
                error('x values must be positive for semilog plot.');
            end
            y = slope * (log10(x) - mean(log10(x))) + y_center;
        case 'loglog'
             % Ensure x values are positive for log10
            if any(x <= 0)
                error('x values must be positive for loglog plot.');
            end
            y = 10.^(slope * (log10(x) - mean(log10(x)))) * y_center;
    end

    plot(x, y, 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); % Exclude from legend
    hold on; % Ensure hold is still on

    % Determine text position and alignment
    if ~isempty(options.TextLocation)
        % Use specified coordinates
        text_x = options.TextLocation(1);
        text_y = options.TextLocation(2);
        % Use centered alignment as a reasonable default for specific coordinates
        horz_align = 'center';
        vert_align = 'middle';
    else
        % Default to middleBottom placement if no coordinates are given
        % Calculate y_offset for default placement
        if max(y) == min(y)
            y_offset = 0.05 * abs(mean(y));
            if y_offset == 0, y_offset = 0.2; end
        else
            y_range = max(y) - min(y);
            y_offset = 0.4 * y_range; % Adjust this factor as needed
        end
        text_x = mean(x);
        text_y = mean(y) - y_offset; % Position below the middle
        horz_align = 'center';
        vert_align = 'top'; % Align top edge of text with text_y
    end

    % Place the text label
    slope_string = rats(slope); % Convert slope to fractional string
    text(text_x, text_y, ['$' slope_string '$'], ...
        'HorizontalAlignment', horz_align, 'VerticalAlignment', vert_align, ...
        'FontSize', 14, 'Interpreter', 'latex');

end