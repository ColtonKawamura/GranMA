function plotEigenmodeDamped(x0, y0, eigen_vectors, mode_to_plot, options)

    arguments
        x0 (:,1) double
        y0 (:,1) double
        eigen_vectors (:,:) double
        mode_to_plot (1,1) double
        options.undamped (1,1) logical = false
    end


    if mode_to_plot > size(eigen_vectors, 2) || mode_to_plot < 1
        error('Invalid mode_to_plot. It must be a valid column index in eigen_vectors.');
    end
     
    N = length(x0); % number of particles

    % grab the egienmode(vector) each column is a mode (eigen_vector = 2N by 4N because each particle has xy coordiantes (rows) and each particle has a  complex conjugate imaginary + and - damped component (columns))
    mode_vector = eigen_vectors(:, mode_to_plot);

    % grab the real, imaginary, and absolute components of the mode
    mode_real = real(mode_vector);
    mode_imag = imag(mode_vector);
    mode_abs = abs(mode_vector);

    %get the dispmment components. Assume odd indices = x, even indices = y)
    dx_real = mode_real(1:2:end); % starts at 1, goes to end, increments by 2 aka only grabs the odd indices
    dy_real = mode_real(2:2:end); % starts at 2, goes to end, increments by 2 aka only grabs the even indices
    dx_imag = mode_imag(1:2:end);
    dy_imag = mode_imag(2:2:end);
    dx_abs = mode_abs(1:2:end);
    dy_abs = mode_abs(2:2:end);
    
    % normalize so it's no so big and wierd looking
    % max_displacement = max([dx_abs; dy_abs]);
    % if max_displacement > 0
    %     dx_real = dx_real / max_displacement;
    %     dy_real = dy_real / max_displacement;
        
    %     dx_imag = dx_imag / max_displacement;
    %     dy_imag = dy_imag / max_displacement;
        
    %     dx_abs = dx_abs / max_displacement;
    %     dy_abs = dy_abs / max_displacement;
    % end

    %  This is the combination of all the stuff that commmented out below:

    % Combined plot
    figure(1); clf;
    hold on;
    quiver(x0, y0, dx_real, dy_real, 0.5, 'b', 'LineWidth', 1.5); % Blue for real part
    quiver(x0, y0, dx_imag, dy_imag, 0.5, 'r', 'LineWidth', 1.5); % Red for imaginary part
    quiver(x0, y0, dx_abs, dy_abs, 0.5, 'g', 'LineWidth', 1.5);   % Green for magnitude
    % scatter(x0, y0, 100 * mode_abs(1:2:end), 'g', 'filled')

    % plot(x0, y0, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');

    title(['\textbf{Eigenmode } $', num2str(mode_to_plot), '$'], 'Interpreter', 'latex');
    xlabel('$x$ position', 'Interpreter', 'latex');
    ylabel('$y$ position', 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    
    legend({'Real', 'Imaginary', 'Magnitude'},'Interpreter', 'latex', 'Location', 'best');
    
    axis equal;
    grid on;
end
