function [matching_omega_gamma_list, loop_mean_wavespeed_list] = plotWavespeed(data, gamma_value, mean_diameter, plot_flag)
    % PLOTWAVESPEED Creates wavespeed plots for granular media
    % 
    % Inputs:
    %   data - Structure containing simulation data loaded from converted MAT file
    %   gamma_value - Damping value to filter by
    %   mean_diameter - Optional mean diameter for normalization (default: 1.0)
    %   plot_flag - Boolean to control plotting (default: true)
    %
    % Outputs:
    %   matching_omega_gamma_list - List of omega*gamma values
    %   loop_mean_wavespeed_list - List of wavespeed values
    
    % Handle optional parameters
    if nargin < 3
        mean_diameter = 1.0;
    end
    if nargin < 4
        plot_flag = true;
    end
    
    % Initialize outputs
    matching_omega_gamma_list = [];
    loop_mean_wavespeed_list = [];

    % Find index closest to desired gamma value
    [~, closest_gamma_index] = min(abs(data.gamma - gamma_value));
    closest_gamma_value = data.gamma(closest_gamma_index);
    
    % Find indices matching the closest gamma value
    matching_gamma_indices = find(data.gamma == closest_gamma_value);
    plot_gamma = gamma_value;

    % Get list of unique pressures for the matching gamma data
    pressure_list = unique(data.pressure(matching_gamma_indices));
    pressure_list = sort(pressure_list);
    plot_pressure = pressure_list;

    % Define the plot limits to match the 1D theory plot curves
    theory_x = (3E-4:1E-5:3);
    theory_y = 1 ./ sqrt(2) .* ((1 + theory_x.^2) .* (1 + sqrt(1 + theory_x.^2))).^(0.5) ./ (1 + theory_x.^2);

    if plot_flag
        % Initialize the figure
        figure_wavespeed = figure;
        loglog(theory_x, 1./theory_y, 'k', 'DisplayName', '1-D Theory');
        hold on;
        xlabel('$\hat{\omega}\hat{\gamma}$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$\hat{c}$', 'FontSize', 20, 'Interpreter', 'latex');
        set(gca, 'XScale', 'log');
        set(get(gca, 'ylabel'), 'rotation', 0);
        grid on;
        box on;
    end

    % Normalize the pressure values for color mapping
    log_pressure = log(pressure_list);
    normalized_variable = (log_pressure - min(log_pressure)) / (max(log_pressure) - min(log_pressure));

    % Create a line for each pressure value
    for p = 1:length(pressure_list)
        pressure_value = pressure_list(p);
        
        % Assign a color based on normalized pressure
        marker_color = [normalized_variable(p), 0, 1-normalized_variable(p)];
        
        % Find indices matching both gamma and current pressure
        matching_pressure_indices = find(data.gamma == closest_gamma_value & ...
                                        data.pressure == pressure_value);
        
        % Get unique omega_gamma values for current pressure
        omega_gamma_values = data.omega_gamma(matching_pressure_indices);
        unique_omega_gamma = unique(omega_gamma_values);
        unique_omega_gamma = sort(unique_omega_gamma);
        
        % Initialize wavespeed list for this pressure
        loop_mean_wavespeed_list = zeros(size(unique_omega_gamma));
        
        % For each omega_gamma value, compute mean wavespeed
        for og = 1:length(unique_omega_gamma)
            omega_gamma_value = unique_omega_gamma(og);
            
            % Find indices matching gamma, pressure, and current omega_gamma
            matching_og_indices = find(data.gamma == closest_gamma_value & ...
                                      data.pressure == pressure_value & ...
                                      data.omega_gamma == omega_gamma_value);
            
            % Compute mean wavespeed for all matching seeds
            wavespeed_values = zeros(length(matching_og_indices), 1);
            for i = 1:length(matching_og_indices)
                idx = matching_og_indices(i);
                % Check if wavespeed_x is a cell array or a regular array
                if iscell(data.wavespeed_x)
                    wavespeed_values(i) = data.wavespeed_x{idx};
                else
                    wavespeed_values(i) = data.wavespeed_x(idx);
                end
            end
            
            % Store the mean value
            loop_mean_wavespeed_list(og) = mean_diameter * mean(wavespeed_values);
        end
        
        % Filter data to include only points where omega_gamma <= 2*gamma_value
        valid_indices = unique_omega_gamma <= gamma_value*2;
        matching_omega_gamma_list = unique_omega_gamma(valid_indices);
        loop_mean_wavespeed_list = loop_mean_wavespeed_list(valid_indices);
        
        % Create legend label
        legend_label = sprintf('$ \\hat{P} = %.3f, \\hat{\\gamma} = %.3f$', ...
                              pressure_value, closest_gamma_value);
        
        % Plot if requested
        if plot_flag
            figure(figure_wavespeed);
            set(gca, 'YScale', 'log');
            plot(matching_omega_gamma_list, loop_mean_wavespeed_list, '-o', ...
                'MarkerFaceColor', marker_color, 'Color', marker_color, ...
                'DisplayName', legend_label);
        end
    end
    
    % Add legend
    if plot_flag
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    end
end