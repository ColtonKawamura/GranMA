export
    plotGausWavespeed2d,
    plotGausWavenumber,
    plotGausAttenuation2d


function plotGausWavespeed2d(simulation_data; plot=true)  
    loop_mean_wavespeed_list = []
    dt = 0.0157

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        % loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        loglog($(theory_x), $(theory_x).^.5, 'k:', 'DisplayName', 'Slope 1/2');
        xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\hat{c} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each omega value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, simulation_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_wavespeed_list = Float64[];
        loop_std_wavespeed_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))

        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            loop_mean_wavespeed = mean(filter(x -> x > 0, [entry.wavespeed for entry in matching_omega_data])) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            loop_std_wavespeed =  std(entry.wavespeed for entry in matching_omega_data) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            push!(loop_std_wavespeed_list, loop_std_wavespeed)

            # Append values
            push!(loop_mean_wavespeed_list, loop_mean_wavespeed)
        end


        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.4f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            x = $(matching_omega_list);
            y = $(loop_mean_wavespeed_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(x, y, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            % errorbar(x, y, $(loop_std_wavespeed_list), '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return loop_mean_attenuation_list
end


function plotGausWavenumber(simulation_data; plot=true)
    loop_mean_wavenumber_list = []
    dt = 0.0157

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        % loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        loglog($(theory_x), $(theory_x).^.5, 'k:', 'DisplayName', 'Slope 1/2');
        xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\hat{k} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each omega value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, simulation_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_wavenumber_list = Float64[];
        loop_std_wavenumber_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))

        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            loop_mean_wavespeed = mean(filter(x -> x > 0, [entry.wavespeed for entry in matching_omega_data])) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter .* dt)
            loop_std_wavespeed =  std(entry.wavespeed for entry in matching_omega_data) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            
            loop_mean_wavenumber = omega_value ./ loop_mean_wavespeed 
            loop_std_wavenumber= omega_value ./ loop_std_wavespeed

            push!(loop_std_wavenumber_list, loop_std_wavenumber)

            # Append values
            push!(loop_mean_wavenumber_list, loop_mean_wavenumber)
        end


        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.4f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            x = $(matching_omega_list);
            y = $(loop_mean_wavenumber_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(x, y, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            % errorbar(x, y, $(loop_std_wavenumber_list), '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
end

function plotGausAttenuation2d(simulation_data; plot=true)  
    loop_mean_attenuation_list = []


    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        % loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        loglog($(theory_x), $(theory_x).^.5, 'k:', 'DisplayName', 'Slope 1/2');
        xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each omega value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, simulation_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];
        loop_std_attenuation_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))

        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean(filter(x -> x > 0, [entry.attenuation for entry in matching_omega_data])) ./ omega_value
            loop_std_alphaoveromega =  std(entry.attenuation for entry in matching_omega_data) ./ omega_value
            push!(loop_std_attenuation_list, loop_std_alphaoveromega)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
        end


        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.3f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(matching_omega_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(omega_gamma, mean_attenuation_x, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label)
            % errorbar(omega_gamma, mean_attenuation_x, $(loop_std_attenuation_list), '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return loop_mean_attenuation_list
end