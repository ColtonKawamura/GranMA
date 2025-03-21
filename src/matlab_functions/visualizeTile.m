function visualizeTile(packing_name)
    % visualizeTile('in/2d_poly_20by20/2D_N12000_P0.001_Width20_Seed1.mat)
        load(packing_name);
    x_mult = 3;
    y_mult = 3;
    N_original = N;
    Dn_original = Dn;
    x_original = x;
    y_original = y;

    N_repeated = x_mult;
    x_repeated = [];
    y_repeated = [];
    Dn_repeated = [];

    for i = 0:N_repeated-1
        x_shifted = x + i * Lx;
        y_shifted = y + i * Ly;
        x_repeated = [x_repeated, x_shifted];
        y_repeated = [y_repeated, y];
        Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
    end

    N = N * N_repeated;
    x = x_repeated;
    y = y_repeated;
    Dn = Dn_repeated; % Set Dn to the repeated diameters

    % Time for y
    N_repeated = y_mult;
    x_repeated = [];
    y_repeated = [];
    Dn_repeated = [];

    for i = 0:N_repeated-1
        y_shifted = y + i * Ly;
        x_repeated = [x_repeated, x];
        y_repeated = [y_repeated, y_shifted];
        Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
    end

    N = N * N_repeated;
    x = x_repeated;
    y = y_repeated;
    Dn = Dn_repeated; % Set Dn to the repeated diameters
    % figure;
    % hold on;
    % axis equal;
    % for np = 1:N
    %     rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], 'Curvature', [1, 1], 'EdgeColor', 'b');
    % end
    % axis([0, N_repeated * Lx, 0, Ly]);
    % hold off;
    % pause
    figure;
    hold on;
    axis equal;
    axis([Lx, 2* Lx, Ly, Ly*2]);

    % Loop through each particle and plot its rectangle
    for np = 1:N
        % Compute the position for the rectangle using the center coordinates (x, y)
        % and the diameter Dn (width and height).
        rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], ...
            'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5); % Optional LineWidth for better visibility
        

        % Wrap the "cell" index around the total number of particles
        label = mod(np-1, N_original) + 1;  % This will cycle labels from 1 to N
        text(x(np), y(np), num2str(label), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
    end
    % Update the figure display
    drawnow;
    hold off;
end
