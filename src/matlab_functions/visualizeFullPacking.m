function visualizeFullPacking(packing_name, x_limit, y_limit)
% visualizeFullPacking('in/2d_poly_20by20/2D_N12000_P0.001_Width20_Seed1.mat', 20, 20)
    load(packing_name);
    figure;
    hold on;
    axis equal;
    axis([0, x_limit, 0, y_limit]);

    % Loop through each particle and plot its rectangle
    for np = 1:N
        % Compute the position for the rectangle using the center coordinates (x, y)
        % and the diameter Dn (width and height).
        rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], ...
            'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5); % Optional LineWidth for better visibility
    end
    % Update the figure display
    drawnow;
    hold off;
end