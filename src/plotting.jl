export
    plotGausWavespeed2d,
    plotGausWavenumber2d,
    plotGausAttenuation2d,
    plotGausAttenuationK2d,
    plot_ωγ_attenuation_2d,
    plot_ωγ_wavespeed_2d,
    plotGuess,
    plotAmp


function plotGausWavespeed2d(simulation_data; plot=true)  
    loop_mean_wavespeed_list = []
    pressure_out = Float64[]
    wavespeed_out = Float64[]

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
            loop_mean_wavespeed = mean(filter(x -> x > 0, [entry.wavespeed for entry in matching_omega_data]))
            loop_std_wavespeed =  std(entry.wavespeed for entry in matching_omega_data) 
            # loop_mean_wavespeed = mean(filter(x -> x > 0, [entry.wavespeed for entry in matching_omega_data])) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            # loop_std_wavespeed =  std(entry.wavespeed for entry in matching_omega_data) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            push!(loop_std_wavespeed_list, loop_std_wavespeed)

            # Append values
            push!(loop_mean_wavespeed_list, loop_mean_wavespeed)
        end
        mean_wavespeed = mean(loop_mean_wavespeed_list)
        std_wavespeed = std(loop_mean_wavespeed_list)
        lower_bound = mean_wavespeed - 2 * std_wavespeed
        upper_bound = mean_wavespeed + 2 * std_wavespeed
        cleaned_wavespeed = filter(x -> x >= lower_bound && x <= upper_bound, loop_mean_wavespeed_list)
        mean_cleaned_wavespeed = mean(cleaned_wavespeed)
        println("$pressure_value")
        pressure_out = push!(pressure_out, pressure_value)
        wavespeed_out = push!(wavespeed_out, mean_cleaned_wavespeed)

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
    return pressure_out, wavespeed_out
end

#### Gaus Plots

function plotGausWavenumber2d(simulation_data; plot=true)
    loop_mean_wavenumber_list = []
    wavenumber_out = Float64[]
    pressure_out, wavespeed_out = plotGausWavespeed2d(simulation_data; plot=false) 

    # Get a list of unique input pressures
    # pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    pressure_list = sort(unique(pressure for pressure in pressure_out))
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
        wavespeed_at_pressure = wavespeed_out[idx]

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
            # loop_mean_wavespeed = mean(filter(x -> x > 0, [entry.wavespeed for entry in matching_omega_data])) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter)
            loop_std_wavespeed =  std(entry.wavespeed for entry in matching_omega_data) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            loop_mean_wavespeed = wavespeed_at_pressure

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
    return wavenumber_out, pressure_out
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
        loglog($(theory_x), $(theory_x).^.5 ./ 15, 'k:', 'DisplayName', 'Slope 1/2');
        loglog($(theory_x), $(theory_x).^1.5 ./ 15, 'k', 'DisplayName', 'Slope 3/2');
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
            plot(omega_gamma, mean_attenuation_x, '-o', 'Color', marker_color, 'DisplayName', legend_label)
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

function plotGausAttenuationK2d(simulation_data; plot=true)  
    loop_mean_attenuation_list = []
    pressure_out, wavespeed_out = plotGausWavespeed2d(simulation_data; plot=false)

    # Get a list of unique input pressures
    # pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    pressure_list = sort(unique(pressure for pressure in pressure_out))
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
        loglog($(theory_x), $(theory_x).^.5 ./ 15, 'k:', 'DisplayName', 'Slope 1/2');
        loglog($(theory_x), $(theory_x).^1.5 ./ 15, 'k', 'DisplayName', 'Slope 3/2');
        xlabel('\$\\hat{k}\$', "FontSize", 20, "Interpreter", "latex");
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
        wavenumber_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        # matching_wavenumber_list = sort(unique([entry.wavenumber for entry in matching_pressure_data]))
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))
        # matching_wavenumber_list = wavenumber_out
        
        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            # matching_wavenumber_data = filter(entry -> entry.wavenumber == wavenumber_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            # omega_value = matching_wavenumber_data[1].omega
            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean(filter(x -> x > 0, [entry.attenuation for entry in matching_omega_data]))  ./ omega_value
            loop_std_alphaoveromega =  std(entry.attenuation for entry in matching_omega_data) ./ omega_value
            push!(loop_std_attenuation_list, loop_std_alphaoveromega)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
            loop_wavenumber_ = omega_value / wavespeed_out[idx]
            push!(wavenumber_list, loop_wavenumber_)
        end
        
        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.4f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(wavenumber_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(omega_gamma, mean_attenuation_x, '-|', 'Color', marker_color, 'DisplayName', legend_label)
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

#### Regular Plots

function plot_ωγ_attenuation_2d(simulation_data, gamma_value, mean_diameter; plot=true)
    # Initialize outputs
    matching_omega_gamma_list = []
    loop_mean_attenuation_list = []

    # filter the data based on those that are close to gamma_value
    closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
    closest_gamma_value = simulation_data[closest_gamma_index].gamma
    matching_gamma_data = filter(entry -> entry.gamma == closest_gamma_value, simulation_data)
    plot_gamma = gamma_value

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_gamma_data])) # goes through each entry of simulation_data and get the P value at that entry
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
        loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        hold on;
        xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
        set(gca, 'XScale', 'log');
        set(get(gca, 'ylabel'), 'rotation', 0);
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each gamma value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))

        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current pressure value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean_diameter .* mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]

        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.3f, \\hat{\\gamma} = %.3f\$", pressure_value, gamma_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(matching_omega_gamma_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            plot_gamma = $(plot_gamma);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(omega_gamma, mean_attenuation_x, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
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
    return matching_omega_gamma_list, loop_mean_attenuation_list
end

function plot_ωγ_wavespeed_2d(simulation_data, gamma_value) # Need to fix lgend

    # filter the data based on those that are close to gamma_value
    closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
    closest_gamma_value = simulation_data[closest_gamma_index].gamma
    matching_gamma_data = filter(entry -> entry.gamma == closest_gamma_value, simulation_data)
    plot_gamma = gamma_value

    # Get a list of unique input pressures
    pressure_list = unique([entry.pressure for entry in matching_gamma_data]) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = 1 ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(0.5) ./ (1 .+ theory_x.^2);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    # Intialized the plots to loop over
    mat"""
    figure_wavespeed = figure;
    loglog($(theory_x), 1./$(theory_y), 'k', 'DisplayName', '1-D Theory'), hold on
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$\\hat{c} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    """

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each gamma value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_wavespeed_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))

        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current pressure value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            @bp
            # Get the mean over all seeds
            loop_mean_wavespeed = mean(entry.:wavespeed_x for entry in matching_omega_gamma_data)
            @bp
            # Append values
            push!(loop_mean_wavespeed_list, loop_mean_wavespeed)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        loop_mean_wavespeed_list = loop_mean_wavespeed_list[valid_indices]
        
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        mean_wavespeed_x = $(loop_mean_wavespeed_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = sprintf('Wavespeed X = %.2f, Gamma = %.2f (Aspect Ratio)', $(pressure_value), $(plot_gamma));
        % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
        pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Wavespeed X)', $(pressure_value), $(plot_gamma));
        
        figure(figure_wavespeed);
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_wavespeed_x, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
        """
    end

    # Add legends to the plots
    mat"""
    % Add legends to the MATLAB plots
    figure(figure_wavespeed);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end

function plotGuess(simulation_data; plot=true)  # This is expirimental
    loop_mean_attenuation_list = []
    pressure_out, wavespeed_out = plotGausWavespeed2d(simulation_data; plot=false)

    # Get a list of unique input pressures
    # pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    pressure_list = sort(unique(pressure for pressure in pressure_out))
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
        %loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        %loglog($(theory_x), $(theory_x).^.5 ./ 15, 'k:', 'DisplayName', 'Slope 1/2');
        %loglog($(theory_x), $(theory_x).^1.5 ./ 15, 'k', 'DisplayName', 'Slope 3/2');
        xlabel('\$\\hat{k}\\hat{P}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$\\hat{\\alpha} \$', "FontSize", 20, "Interpreter", "latex");
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
        wavenumber_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        # matching_wavenumber_list = sort(unique([entry.wavenumber for entry in matching_pressure_data]))
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))
        # matching_wavenumber_list = wavenumber_out
        
        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            # matching_wavenumber_data = filter(entry -> entry.wavenumber == wavenumber_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            # omega_value = matching_wavenumber_data[1].omega
            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean(filter(x -> x > 0, [entry.attenuation for entry in matching_omega_data]))  #./ omega_value
            loop_std_alphaoveromega =  std(entry.attenuation for entry in matching_omega_data) #./ omega_value
            push!(loop_std_attenuation_list, loop_std_alphaoveromega)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
            loop_wavenumber_ = omega_value / wavespeed_out[idx] #* pressure_value
            push!(wavenumber_list, loop_wavenumber_)
        end
        
        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.4f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(wavenumber_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(omega_gamma, mean_attenuation_x, '-|', 'Color', marker_color, 'DisplayName', legend_label)
            xline($(pressure_value) * .5)
            yline(.12)
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

# Single Simulaiton plots

function plotAmp(filtered_data)
    x = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    y = filtered_data[1].amplitude_vector_x
    display_name = @sprintf("\$ A_{||}(x) \$")
    mat"""
    figure
    scatter($(x), $(y), "DisplayName", $(display_name))
    set(gca, 'YScale', 'log')
    grid on
    xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
    ylabel("\$ A(x) \$", "Interpreter", 'latex', "FontSize", 15)
    set(get(gca, 'ylabel'), 'rotation', 0);
    box on
    hold on 
    """
    x = filtered_data[1].initial_distance_from_oscillation_output_y_fft
    y = filtered_data[1].amplitude_vector_y
    display_name = @sprintf("\$ A_\\perp(x) \$")
    mat"""
    scatter($(x), $(y), "DisplayName", $(display_name))
    set(gca, 'YScale', 'log')
    grid on
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end