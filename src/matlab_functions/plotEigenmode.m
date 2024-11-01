function plotEigenmode(x0, y0, eigen_vectors, mode_to_plot)
    if mode_to_plot > size(eigen_vectors, 2) || mode_to_plot < 1
        error('Invalid mode_to_plot. It must be a valid column index in eigen_vectors.');
    end
    
    N = length(x0);

    mode_positions = eigen_vectors(:, mode_to_plot); % grab the mode (column of positions)
    
    dx = zeros(N, 1);
    dy = zeros(N, 1);
    
    % Populate dx and dy for each particle based on eigen_vectors
    for i = 1:N
        dx(i) = mode_positions(2*i-1) - x0(i); % x displacement
        dy(i) = mode_positions(2*i) - y0(i);   % y displacement
    end
    
    quiver(x0, y0, dx, dy, 0, 'AutoScale', 'off');
    title(['Quiver plot for eigenmode ', num2str(mode_to_plot)]);
    xlabel('X Position');
    ylabel('Y Position');
    axis equal;
    grid on;
end
