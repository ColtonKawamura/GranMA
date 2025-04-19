function [mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase] = plotMeanField(filtered_data, options)
    % Rewritten from Julia to MATLAB
    arguments
        filtered_data (1,1) struct % Required argument: scalar struct containing simulation data fields
        options.plot (1,1) logical = true % Optional flag to generate plots
        options.shear (1,1) logical = false % Optional flag for shear vs compression
        options.z (1,1) logical = false % plots z direction instead
    end

    % Access fields directly from the scalar struct
    % Extract data from cell arrays using {1} assuming filterData reduced to one entry
    if options.shear
        amplitude_vector = filtered_data.amplitude_vector_y{1}; % Extract from cell
        attenuation = filtered_data.alphaoveromega_y; % Numeric field
        distance_from_wall = filtered_data.initial_distance_from_oscillation_output_y_fft{1}; % Extract from cell
        omega = filtered_data.omega; % Numeric field
        legend_para_prime = '$ \hat{y}'' $'; % LaTeX format for legend
        legend_para = '$ \hat{y} $';
        legend_para_mean = '$ \hat{\overline{y}} $'; % Double backslash for overline
        legend_perp = '$ \hat{x} $';
    else
        amplitude_vector = filtered_data.amplitude_vector_x{1}; % Extract from cell
        attenuation = filtered_data.alphaoveromega_x; % Numeric field
        distance_from_wall = filtered_data.initial_distance_from_oscillation_output_x_fft{1}; % Extract from cell
        omega = filtered_data.omega; % Numeric field
        legend_para_prime = '$ \hat{x}'' $';
        legend_para = '$ \hat{x} $';
        legend_para_mean = '$ \hat{\overline{x}} $';
        if options.z
            legend_perp = '$ \hat{z} $';
        else
            legend_perp = '$ \hat{y} $';
        end
    end

    % Ensure vectors are column vectors for polyfit
    distance_from_wall = distance_from_wall(:);
    amplitude_vector = amplitude_vector(:);

    % Calculate mean field amplitude based on initial pressure and attenuation
    A = filtered_data.pressure / 100;
    % mean_field_amp = A * exp(-attenuation * omega * distance_from_wall); % Original calculation

    % Fit amplitude data to get attenuation (alternative approach)
    % Use polyfit on log scale for exponential decay
    % Now abs() should work on the numeric vector extracted from the cell
    coefficients_amp = polyfit(distance_from_wall, log(abs(amplitude_vector)), 1);
    fitted_attenuation = coefficients_amp(1); % Slope corresponds to -attenuation*omega (or similar)
    intercept_attenuation = coefficients_amp(2); % Intercept related to initial amplitude log(A)
    mean_field_amp_fitted = exp(intercept_attenuation) * exp(fitted_attenuation * distance_from_wall);

    % Calculate prime field amplitude (fluctuations around the mean)
    if options.shear
        y = filtered_data.amplitude_vector_y{1}(:); % Extract & ensure column
        x_parra = filtered_data.initial_distance_from_oscillation_output_y_fft{1}(:); % Extract & ensure column
        prime_field_amp = abs(y - mean_field_amp_fitted);
    else
        y = filtered_data.amplitude_vector_x{1}(:); % Extract & ensure column
        x_parra = filtered_data.initial_distance_from_oscillation_output_x_fft{1}(:); % Extract & ensure column
        prime_field_amp = abs(y - mean_field_amp_fitted);
    end

    driving_amp = filtered_data.pressure / 100;

    % Get perpendicular component data
    if options.shear
        x_perp = filtered_data.initial_distance_from_oscillation_output_x_fft{1}(:); % Extract & ensure column
        y_perp = filtered_data.amplitude_vector_x{1}(:); % Extract & ensure column
    else
        if options.z
            x_perp = filtered_data.initial_distance_from_oscillation_output_y_fft{1}(:); % Extract & ensure column
            y_perp = filtered_data.amplitude_vector_y{1}(:); % Extract & ensure column
        else
            x_perp = filtered_data.initial_distance_from_oscillation_output_z_fft{1}(:); % Extract & ensure column
            y_perp = filtered_data.amplitude_vector_z{1}(:); % Extract & ensure column
        end
    end

    % --- Amplitude Plotting ---
    if options.plot
        figure; % Create a new figure for amplitude plots
        hold on; % Hold on from the start
        plotHandles = []; % Initialize empty array for handles
        plotLabels = {}; % Initialize empty cell array for labels

        % Plot Prime Field Amplitude first (Explicit Blue)
        if ~isempty(x_parra) && ~isempty(prime_field_amp)
            h_prime = scatter(x_parra, prime_field_amp / driving_amp, '*', 'MarkerEdgeColor', [0 0.4470 0.7410]); % Default Blue
            if ~isempty(h_prime)
                plotHandles(end+1) = h_prime(1); % Take only the first handle
                plotLabels{end+1} = legend_para_prime;
            end
        end

        % Plot Total Amplitude second (Explicit Orange)
        if ~isempty(x_parra) && ~isempty(y)
            h_total = scatter(x_parra, abs(y) / driving_amp, 'o', 'MarkerEdgeColor', [0.8500 0.3250 0.0980]); % Default Orange
            if ~isempty(h_total)
                plotHandles(end+1) = h_total(1); % Take only the first handle
                plotLabels{end+1} = legend_para;
            end
        end

        % Plot Perpendicular Amplitude third (Explicit Purple)
        if ~isempty(x_perp) && ~isempty(y_perp)
            h_perp = scatter(x_perp, abs(y_perp) / driving_amp, '+', 'MarkerEdgeColor', [0.4940 0.1840 0.5560]); % Default Purple
            if ~isempty(h_perp)
                plotHandles(end+1) = h_perp(1); % Take only the first handle
                plotLabels{end+1} = legend_perp;
            end
        end

        % Plot Fitted Mean Field Amplitude last (Explicitly black)
        if ~isempty(x_parra) && ~isempty(mean_field_amp_fitted)
            h_mean = plot(x_parra, mean_field_amp_fitted / driving_amp, ':', 'Color', 'k', 'LineWidth', 3); % Keep black dotted line
            if ~isempty(h_mean)
                plotHandles(end+1) = h_mean(1); % Take only the first handle
                plotLabels{end+1} = legend_para_mean;
            end
        end

        set(gca, 'YScale', 'log');
        grid on;
        xlabel('$ x_0 $', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$A(x_0)$', 'Interpreter', 'latex', 'FontSize', 15);
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on;
        title('Amplitude Components');

        % Create legend using collected handles and labels if any plots were made
        if ~isempty(plotHandles)
            legend(plotHandles, plotLabels, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);
        end
        set(gca, 'FontSize', 15);
        hold off;
    end

    % --- Phase Calculation ---
    if options.shear
        phase = filtered_data.unwrapped_phase_vector_y{1}(:); % Extract & ensure column
        distance_from_wall_phase = filtered_data.initial_distance_from_oscillation_output_y_fft{1}(:); % Extract & ensure column
        wavenumber = filtered_data.wavenumber_y;
    else
        phase = filtered_data.unwrapped_phase_vector_x{1}(:); % Extract & ensure column
        distance_from_wall_phase = filtered_data.initial_distance_from_oscillation_output_x_fft{1}(:); % Extract & ensure column
        wavenumber = filtered_data.wavenumber_x;
    end
    phase_mod = mod(phase, 2*pi); % Wrap phase to [0, 2*pi)

    % Unwrap scattered phase data and fit linear trend for mean field phase
    % Modified call to handle potential duplicate distances
    [distance_from_wall_unique_sorted, unwrapped_phase_mean_sorted] = unwrapScattered(distance_from_wall_phase, phase); % Using original phase before mod

    % Fit a line to the unwrapped mean phase vs unique distance
    fitline_coeffs = polyfit(distance_from_wall_unique_sorted, unwrapped_phase_mean_sorted, 1);
    % Evaluate the fitted line at the unique sorted distances
    mean_field_phase_unwrapped_trend = polyval(fitline_coeffs, distance_from_wall_unique_sorted);

    % Calculate prime field phase (fluctuations around the mean trend)
    % Interpolate the mean field trend back to the original distances
    mean_field_phase_interp = interp1(distance_from_wall_unique_sorted, mean_field_phase_unwrapped_trend, distance_from_wall_phase, 'linear', 'extrap');
    prime_field_phase = abs(phase - mean_field_phase_interp); % Difference from unwrapped original phase
    prime_field_phase_mod = mod(prime_field_phase, 2*pi); % Wrap prime phase

    % Wrap the mean phase trend for plotting
    mean_field_phase = mod(mean_field_phase_unwrapped_trend, 2*pi);

    % Get perpendicular phase component data
    if options.shear
        x_perp_phase = filtered_data.initial_distance_from_oscillation_output_x_fft{1}(:); % Extract & ensure column
        y_perp_phase = filtered_data.unwrapped_phase_vector_x{1}(:); % Extract & ensure column
    else
        if options.z
            x_perp_phase = filtered_data.initial_distance_from_oscillation_output_z_fft{1}(:); % Extract & ensure column
            y_perp_phase = filtered_data.unwrapped_phase_vector_z{1}(:); % Extract & ensure column
        else
            x_perp_phase = filtered_data.initial_distance_from_oscillation_output_y_fft{1}(:); % Extract & ensure column
            y_perp_phase = filtered_data.unwrapped_phase_vector_y{1}(:); % Extract & ensure column
        end
    end
    y_perp_phase_mod = mod(y_perp_phase, 2*pi);

    % --- Phase Plotting ---
    if options.plot
        figure; % Create a new figure for phase plots

        % Plot Prime Field Phase
        scatter(distance_from_wall_phase, prime_field_phase_mod, '*', 'DisplayName', '$\phi_{||}''$'); % Use prime_field_phase_mod
        hold on;
        grid on;
        xlabel('$ x $', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$ \phi(x) $', 'Interpreter', 'latex', 'FontSize', 15);
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on;
        title('Phase Components'); % Add a title

        % Plot Total Phase (wrapped)
        scatter(distance_from_wall_phase, phase_mod, 'o', 'DisplayName', '$\phi_{||}$');

        % Plot Mean Field Phase (wrapped) - Use unique sorted distances
        scatter(distance_from_wall_unique_sorted, mean_field_phase, 'v', 'MarkerEdgeColor', 'k', 'SizeData', 15, 'DisplayName', '$\overline{\phi}_{||}$'); % Smaller markers for mean

        % Plot Perpendicular Phase (wrapped)
        scatter(x_perp_phase, y_perp_phase_mod, '+', 'DisplayName', '$\phi_\perp$');

        legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);
        set(gca, 'FontSize', 15);
        hold off;
    end

    % Assign output variables (using fitted amplitude and calculated phase)
    mean_field_amp = mean_field_amp_fitted;
    % mean_field_phase is already assigned (wrapped version)
    % prime_field_amp is already assigned
    % prime_field_phase is already assigned (unwrapped difference)

end % End of function

% Helper funciton
function [x_unique_sorted, y_unwrapped_mean_sorted] = unwrapScattered(x, y_wrapped)
    % Handle duplicate x values by averaging y values, then unwrap.
    [x_unique, ~, ic] = unique(x); % Find unique x values and indices for grouping

    % Calculate the mean of y_wrapped for each unique x value
    % Ensure y_wrapped is a column vector for accumarray
    y_mean_at_unique_x = accumarray(ic, y_wrapped(:), [], @mean); 

    % Sort the unique x values and corresponding mean y values
    [x_unique_sorted, sort_idx] = sort(x_unique);
    y_mean_sorted_wrapped = y_mean_at_unique_x(sort_idx);

    % Unwrap the sorted mean y values
    y_unwrapped_mean_sorted = unwrap(y_mean_sorted_wrapped);
    % Optional: Add check for large jumps if needed
end