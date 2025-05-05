function plotEigenmode(x0, y0, eigenVectors, modeToPlot, options)
% Note: user chooses modde to plot but the vectors come in cojugate pairs
% I shoudl considder just cleaning up the eigenVectors and values with only the real
% eigen freqs before even passing to this....

    arguments
        x0 (:,1) double
        y0 (:,1) double
        eigenVectors (:,:) double
        modeToPlot (1,1) double
        options.damped (1,1) logical = false
    end

    % grab the egienmode(vector) each column is a mode (eigen_vector = 2N by 4N because each particle has xy coordiantes (rows) and each particle has a  complex conjugate imaginary + and - damped component (columns))
    modeVector = eigenVectors(:, modeToPlot);

    %% Rotate each eigenVector in each mode by it's unit phasor
    if options.damped
        % eigenValueKeep = imag(eigenValues) > 0; % keep the positive eigenvalues
        [~, idx]   = max(abs(modeVector)); % get the largest eigenvector this column
        % Turn largest eigenvector into a unit phasor
        unitPhasor =  modeVector(idx) / abs(modeVector(idx));
        modeVector= modeVector(:) / unitPhasor; % rotate all the eigenvectors by the same amount
    end


    if options.damped
        % grab the real, imaginary, and absolute components of the mode
        mode_real = modeVector;
        % mode_imag = imag(modeVector);
        % mode_abs = abs(modeVector);
    else
        % grab the real, imaginary, and absolute components of the mode
        mode_real = real(modeVector);
    end


    if options.damped
        dx_real = modeVector(1:2:end); % starts at 1, goes to end, increments by 2 aka only grabs the odd indices
        dy_real = modeVector(2:2:end); % starts at 2, goes to end, increments by 2 aka only grabs the even indices% 
        % dx_imag = mode_imag(1:2:end);
        % dy_imag = mode_imag(2:2:end);
        % dx_abs = mode_abs(1:2:end);
        % dy_abs = mode_abs(2:2:end);
    else
        %get the dispmment components. Assume odd indices = x, even indices = y)
        dx_real = mode_real(1:2:end); % starts at 1, goes to end, increments by 2 aka only grabs the odd indices
        dy_real = mode_real(2:2:end); % starts at 2, goes to end, increments by 2 aka only grabs the even indices
    end

    figure(1); clf;
    hold on;
    quiver(x0, y0, dx_real, dy_real, 0.5, 'b', 'LineWidth', 1.5); % Blue for real part
    if options.damped
        % quiver(x0, y0, dx_imag, dy_imag, 0.5, 'r', 'LineWidth', 1.5); % Red for imaginary part
        % quiver(x0, y0, dx_abs, dy_abs, 0.5, 'g', 'LineWidth', 1.5);   % Green for magnitude
    end
    % scatter(x0, y0, 100 * mode_abs(1:2:end), 'g', 'filled')
    % plot(x0, y0, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');

    title(['\textbf{Eigenmode } $', num2str(modeToPlot), '$'], 'Interpreter', 'latex');
    xlabel('$x$ position', 'Interpreter', 'latex');
    ylabel('$y$ position', 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    if options.damped
        % legend({'Real', 'Imaginary', 'Magnitude'},'Interpreter', 'latex', 'Location', 'best');
    else
        legend({'Real'},'Interpreter', 'latex', 'Location', 'best');
    end
    axis equal;
    grid on;
end
