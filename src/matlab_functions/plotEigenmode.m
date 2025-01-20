function plotEigenmode(x0, y0, eigen_vectors, mode_to_plot, damped)
    % if system is damped, set damped = true. Otherwise, set damped = false.

    % clunky matlab way to set default values
    if nargin < 5
        damped = false;
    end

    if mode_to_plot > size(eigen_vectors, 2) || mode_to_plot < 1
        error('Invalid mode_to_plot. It must be a valid column index in eigen_vectors.');
    end
     
    N = length(x0);

    if damped == true
        mode_positions = abs(eigen_vectors(:, mode_to_plot)); 
    else
        mode_positions = eigen_vectors(:, mode_to_plot); % grab the mode (column of positions)
    end
    
    dx = zeros(N, 1);
    dy = zeros(N, 1);
    
    % Populate dx and dy for each particle based on eigen_vectors
    for i = 1:N
        dx(i) = mode_positions(2*i-1);% - x0(i); % Particle  x-position are odd rows 2*i-1 = 1,3,4...
        dy(i) = mode_positions(2*i); % - y0(i);   % y-positions are even rows 2*i = 2,4,6...
    end
   figure(1); clf
    quiver(x0, y0, dx, dy)
    title(['Eigenmode ', num2str(mode_to_plot)]);
    % axis equal;
    grid on;
    % ylim([0,y_limit])
    figure(2); clf
     plot(x0,dx, '.'); title("x0 vs dx")
    hold on
    plot(x0,dy, '.') 
end
