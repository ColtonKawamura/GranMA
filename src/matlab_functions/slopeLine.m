function slopeLine(type, slope, x_bounds, y_center, options)
    % SLOPELINE Plot a line with a given slope and center point, with text label placement.
    %   SLOPELINE(TYPE, SLOPE, X_BOUNDS, Y_CENTER) plots a line with default text placement.
    %   SLOPELINE(..., 'TextPosition', POSITION) specifies the text label position.
    %   POSITION can be 'bottomLeft', 'bottomRight', 'topLeft', 'topRight',
    %   'middleTop', 'middleBottom' (default).
    
    arguments
        type {mustBeMember(type, {'linear', 'semilog', 'loglog'})}
        slope (1,1) double
        x_bounds (1,:) double % Ensure x_bounds is a row vector
        y_center (1,1) double
        options.TextPosition {mustBeMember(options.TextPosition, ...
            {'bottomLeft', 'bottomRight', 'topLeft', 'topRight', 'middleTop', 'middleBottom'})} = 'middleBottom'
    end
    
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
    
    % Handle cases where y might be constant (y_range = 0)
    if max(y) == min(y)
        y_offset = 0.05 * abs(mean(y)); % Use a fraction of the mean y value
        if y_offset == 0 % If mean(y) is also 0, use a small absolute offset
            y_offset = 0.2; 
        end
    else
        y_range = max(y) - min(y);
        y_offset = 0.4 * y_range;
    end
    
    
    % Determine text position and alignment
    text_x = mean(x);
    text_y = mean(y);
    horz_align = 'center';
    vert_align = 'top'; % Default for middleBottom
    
    switch options.TextPosition
        case 'bottomLeft'
            text_x = min(x);
            text_y = y(1) - y_offset;
            horz_align = 'left';
            vert_align = 'top';
        case 'bottomRight'
            text_x = max(x);
            text_y = y(end) - y_offset;
            horz_align = 'right';
            vert_align = 'top';
        case 'topLeft'
            text_x = min(x);
            text_y = y(1) + y_offset;
            horz_align = 'left';
            vert_align = 'bottom';
        case 'topRight'
            text_x = max(x);
            text_y = y(end) + y_offset;
            horz_align = 'right';
            vert_align = 'bottom';
        case 'middleTop'
            text_x = mean(x);
            text_y = mean(y) + y_offset;
            horz_align = 'center';
            vert_align = 'bottom';
        case 'middleBottom' % Default case
            text_x = mean(x);
            text_y = mean(y) - y_offset;
            horz_align = 'center';
            vert_align = 'top';
    end
    
    % Place the text label
    slope_string = rats(slope); % Convert slope to fractional string
    text(text_x, text_y, ['$' slope_string '$'], ...
        'HorizontalAlignment', horz_align, 'VerticalAlignment', vert_align, ...
        'FontSize', 14, 'Interpreter', 'latex');
    
    % Removed the extra end statement if it wasn't needed for the function block