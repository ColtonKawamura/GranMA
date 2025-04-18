function plotAttenuationMult(simulation_data, gamma_values, mean_diameter, varargin)
    % Optional argument parsing (example, 'shear' is not used in the current logic)
    % p = inputParser;
    % addParameter(p, 'shear', false, @islogical);
    % parse(p, varargin{:});
    % shear = p.Results.shear;

    % Handle optional mean_diameter
    if nargin < 3
        mean_diameter = 1.0; % Default mean diameter
    end

    shave = false; % forgot why I put this here, need to fix in the future

    % Define the plot limits to match the 1D theory plot curves
    theory_x = 3E-4:1E-5:3;
    theory_y = theory_x ./ sqrt(2) .* ((1 + theory_x.^2) .* (1 + sqrt(1 + theory_x.^2))).^(-0.5);

    figure;
    loglog(theory_x, theory_y, 'k', 'DisplayName', '1-D Theory');
    hold on;
    % Use single backslashes for LaTeX
    xlabel('$\hat{\omega}\hat{\gamma}$', 'FontSize', 20, 'Interpreter', 'latex'); % Checked: Single backslashes
    ylabel('$ \frac{\hat{\alpha}}{\hat{\omega}} $', 'FontSize', 20, 'Interpreter', 'latex'); % Checked: Single backslashes
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;

    marker_shape_vector = {"-*", "-o", "-v", "-+", "-.", "-x", "-d"};

    % Assuming simulation_data is a struct array with fields: gamma, omega, pressure, omega_gamma, alphaoveromega_x
    all_gammas = [simulation_data.gamma]; % Extract all gamma values once

    for g_idx = 1:length(gamma_values)
        gamma_value_target = gamma_values(g_idx);

        % filter the data based on those that are close to gamma_value_target
        [~, closest_gamma_index] = min(abs(all_gammas - gamma_value_target));
        % Get the actual closest gamma value present in the data
        closest_gamma_value = all_gammas(closest_gamma_index);

        % Filter by the closest actual gamma value
        matching_gamma_indices = find([simulation_data.gamma] == closest_gamma_value);

        % Check if any data matches this gamma
        if isempty(matching_gamma_indices)
            warning('No data found for target gamma %.4f (closest actual: %.4f). Skipping.', gamma_value_target, closest_gamma_value);
            continue;
        end

        % Get the subset of data for this gamma
        % Handle scalar simulation_data case explicitly
        if isscalar(simulation_data)
            % If it's scalar, matching_gamma_indices should be 1 if it matched
            if any(matching_gamma_indices) % Check if the scalar struct's gamma matched
                matching_gamma_data = simulation_data;
            else
                % This case should theoretically be caught by the isempty check above,
                % but handle defensively.
                warning('Scalar simulation_data did not match closest gamma %.4f. Skipping.', closest_gamma_value);
                continue;
            end
        else
            % If it's an array, use the indices found
            matching_gamma_data = simulation_data(matching_gamma_indices);
        end

        % Get a list of unique input pressures for this gamma subset
        pressure_list_full = unique([matching_gamma_data.pressure]);
        pressure_list_full = sort(pressure_list_full); % Ensure sorted

        % Check if any pressures were found
        if isempty(pressure_list_full)
            warning('No pressure data found for gamma %.4f. Skipping.', closest_gamma_value);
            continue;
        end

        % get a range for plotting color from 0 to 1 based on all pressures
        log_pressures = log(pressure_list_full);
        min_log_p = min(log_pressures);
        max_log_p = max(log_pressures);
        if max_log_p == min_log_p || isnan(max_log_p) || isnan(min_log_p) % Handle single pressure or potential NaNs
            normalized_variable = zeros(size(pressure_list_full)); % Avoid division by zero or NaN result
        else
            normalized_variable = (log_pressures - min_log_p) ./ (max_log_p - min_log_p);
        end

        % Create a line for each pressure value found
        for p_idx = 1:length(pressure_list_full) % Iterate through ALL pressures
            pressure_value = pressure_list_full(p_idx);

            % Assign a color based on normalized pressure
            marker_color = [normalized_variable(p_idx), 0, 1-normalized_variable(p_idx)];

            % Only look at data for current pressure value within the gamma-filtered data
            matching_pressure_indices = find([matching_gamma_data.pressure] == pressure_value);

            % Handle scalar matching_gamma_data case within the pressure loop
            if isscalar(matching_gamma_data)
                % If scalar, indices should be 1 if pressure matches, [] otherwise
                if any(matching_pressure_indices) % Check if the scalar struct's pressure matched
                    matching_pressure_data = matching_gamma_data;
                else
                    % This case should be impossible if pressure_list_full was derived correctly
                    warning('Scalar matching_gamma_data pressure mismatch. P=%.4f, gamma=%.4f', pressure_value, closest_gamma_value);
                    continue;
                end
            else
                 % If it's an array, use the indices found
                 matching_pressure_data = matching_gamma_data(matching_pressure_indices);
            end

            if isempty(matching_pressure_data) % This check might be redundant now but keep for safety
                warning('Unexpected: No data found for P=%.4f, gamma=%.4f after initial filtering.', pressure_value, closest_gamma_value);
                continue; % Skip if no data for this specific pressure (shouldn't happen ideally)
            end

            % Initialized vectors for just this pressure
            loop_mean_attenuation_list = [];
            plot_omega_gamma_list = []; % Initialize corresponding x-values

            % Look at a single omega gamma value since each one spans all seeds
            matching_omega_gamma_list = unique([matching_pressure_data.omega_gamma]);
            matching_omega_gamma_list = sort(matching_omega_gamma_list);

            for og_idx = 1:length(matching_omega_gamma_list)
                omega_gamma_value = matching_omega_gamma_list(og_idx);

                % <<< ADD FILTER SIMILAR TO plotAttenuation.m >>>
                if omega_gamma_value > closest_gamma_value * 2
                    continue; % Skip points where omega*gamma is too high
                end

                % Only look at data for current omega_gamma value
                matching_omega_gamma_indices = find([matching_pressure_data.omega_gamma] == omega_gamma_value);

                % Handle scalar matching_pressure_data case within the omega_gamma loop
                if isscalar(matching_pressure_data)
                    % If scalar, indices should be 1 if omega_gamma matches, [] otherwise
                    if any(matching_omega_gamma_indices) % Check if the scalar struct's omega_gamma matched
                        matching_omega_gamma_data = matching_pressure_data;
                    else
                        % This case should be impossible if matching_omega_gamma_list was derived correctly
                        warning('Scalar matching_pressure_data omega_gamma mismatch. P=%.4f, gamma=%.4f, og=%.4f', pressure_value, closest_gamma_value, omega_gamma_value);
                        continue;
                    end
                else
                    % If it's an array, use the indices found
                    matching_omega_gamma_data = matching_pressure_data(matching_omega_gamma_indices);
                end

                % Get the mean over all seeds for this omega_gamma
                % Ensure alphaoveromega_x exists and handle potential missing data if necessary
                if isfield(matching_omega_gamma_data, 'alphaoveromega_x') && ~isempty([matching_omega_gamma_data.alphaoveromega_x])
                     % Check for NaNs within the data for the mean calculation
                     current_alphaoveromega = [matching_omega_gamma_data.alphaoveromega_x];
                     valid_alpha_indices = ~isnan(current_alphaoveromega);
                     if any(valid_alpha_indices) % Only calculate mean if there are valid numbers
                         loop_mean_alphaoveromega = mean_diameter * mean(current_alphaoveromega(valid_alpha_indices));
                         % Append values
                         loop_mean_attenuation_list(end+1) = loop_mean_alphaoveromega;
                         plot_omega_gamma_list(end+1) = omega_gamma_value; % Append corresponding x-value
                     else
                         warning('All alphaoveromega_x values are NaN for P=%.4f, gamma=%.4f, omega_gamma=%.4f. Skipping point.', pressure_value, closest_gamma_value, omega_gamma_value);
                     end
                else
                     warning('Missing or empty alphaoveromega_x field for P=%.4f, gamma=%.4f, omega_gamma=%.4f. Skipping point.', pressure_value, closest_gamma_value, omega_gamma_value);
                end
            end

            % Label uses the *closest* gamma value found in the data
            % Use single backslashes for LaTeX
            pressure_label = sprintf('$ \hat{P}=%.3f, \hat{\gamma}_{actual}=%.3f $', pressure_value, closest_gamma_value); % Checked: Single backslashes

            % Find marker shape based on the *target* gamma value from the input list
            marker_idx = find(gamma_values == gamma_value_target, 1, 'first');
            if isempty(marker_idx) || marker_idx > length(marker_shape_vector)
                marker_shape = '-o'; % Default marker if something goes wrong
                warning('Could not find target gamma %.4f in gamma_values or index out of bounds.', gamma_value_target);
            else
                 marker_shape = marker_shape_vector{marker_idx}; % Use curly braces for cell array access
            end

            % Plotting
            if ~isempty(plot_omega_gamma_list) % Only plot if there is data
                plot( plot_omega_gamma_list, loop_mean_attenuation_list, marker_shape, 'Color', marker_color, 'DisplayName', pressure_label);
            end
        end
    end

    % Add legends to the plots
    % Use single backslashes for LaTeX
    leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    title(leg, '$\hat{P}, \hat{\gamma}_{actual}$') % Checked: Single backslashes, removed extra spaces

end % End of function