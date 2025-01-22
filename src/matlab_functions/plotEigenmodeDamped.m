function plotEigenmodeDamped(x0, y0, eigen_vectors, mode_to_plot)
    % need to have the packing loaded in adddition to this so you can feed it the inital coordinates
    %  They are currently in in/2d_eigen_mode_test
    %  Pick the matchign pressure you picked in eigen_vectors = findInStruct(results, {'pressure', 'damping'}, {.1, .1}, 'eigen_vectors')

    if mode_to_plot > size(eigen_vectors, 2) || mode_to_plot < 1
        error('Invalid mode_to_plot. It must be a valid column index in eigen_vectors.');
    end
     
    N = length(x0);

    mode_positions = abs(eigen_vectors(:, mode_to_plot)); 
    
    dx = zeros(N, 1);
    dy = zeros(N, 1);
    
    % Populate dx and dy for each particle based on eigen_vectors
    for index = 1:N
        dx(index) = mode_positions(2*index-1);% - x0(i); % Particle  x-position are odd rows 2*i-1 = 1,3,4...
        dy(index) = mode_positions(2*index); % - y0(i);   % y-positions are even rows 2*i = 2,4,6...
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
