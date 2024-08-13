include("src/GranMA.jl")

using .GranMA
using CSV
using DataFrames
using Debugger
using MATLAB
using Printf
using Statistics
using MLJ
using MLJModels

# Read the data table
data_frame = CSV.read("out/processed/K100_ellipse_edits.csv", DataFrame)

simulation_data = load_data("out/processed/2d_biredo_K100_W5.jld2")

function plot_ellipse_pdf(ω_value, γ_value; plot=true, simulation_data=simulation_data)

    # Intitialize the outputs of the function
    probabilities_asp = []
    probabilities_rot = []
    plot_bins_asp = []
    plot_bins_rot = []

    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)

    # filter the data based on those that are close to gamma_value
    closest_ω_index = argmin(abs.([idx.omega for idx in matching_γ_data] .- ω_value))
    closest_ω_value = matching_γ_data[closest_ω_index].omega
    matching_ωγ_data = filter(entry -> entry.omega == closest_ω_value, matching_γ_data)

    # Get a list of unique input pressures
    pressure_list = unique([entry.pressure for entry in matching_ωγ_data]) # goes through each entry of simulation_data and get the P value at that entry
    
    # get a range for plotting color from 0 to 1
    if length(pressure_list) != 1
        normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))
    else
        normalized_variable = [1] # protect against divide by zero if length is 1
    end

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_aspect_ratio = figure;
        xlabel('\$ b/a \$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ p(b/a) \$', "FontSize", 20, "Interpreter", "latex");
        set(gca, 'YScale', 'log');
        set(get(gca, 'ylabel'), 'rotation', 0);
        grid on;
        box on;
        hold on

        figure_rotation_angle = figure;
        xlabel('\$ \\theta \$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ p(\\theta) \$', "FontSize", 20, "Interpreter", "latex");
        set(gca, 'YScale', 'log');
        set(get(gca, 'ylabel'), 'rotation', 0);
        grid on;
        box on;
        hold on
        """
    end

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_ωγ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
        
        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # mean(matrix, dims=2) means across rows. Dims = 1 is cross columns, dims =3 is into the screen
        #  If you have an array [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]], splatting converts this into calling hcat([a1, a2, a3], [b1, b2, b3], [c1, c2, c3]).
        mean_asp_rat_counts = mean(hcat([entry.asp_rat_counts for entry in matching_pressure_data]...), dims=2)
        mean_asp_rat_bins = mean(hcat([entry.asp_rat_bins for entry in matching_pressure_data]...), dims=2)
        mean_rot_ang_counts = mean(hcat([entry.rot_ang_counts for entry in matching_pressure_data]...), dims=2)
        mean_rot_ang_bins = mean(hcat([entry.rot_ang_bins for entry in matching_pressure_data]...), dims=2)

        # # Calculate bin widths
        asp_rat_bin_widths = diff(mean_asp_rat_bins, dims=1)
        rot_ang_bin_widths = diff(mean_rot_ang_bins, dims=1)

        # # Normalize counts to get probabilities, considering the bin widths
        total_asp_count = sum(mean_asp_rat_counts .* asp_rat_bin_widths)
        total_rot_count = sum(mean_rot_ang_counts .* rot_ang_bin_widths)
        probabilities_asp = (mean_asp_rat_counts ./ total_asp_count) ./ asp_rat_bin_widths
        probabilities_rot = (mean_rot_ang_counts ./ total_rot_count) ./ rot_ang_bin_widths
        plot_bins_asp = mean_asp_rat_bins[1:end-1]
        plot_bins_rot = mean_rot_ang_bins[1:end-1]
        # probabilities_asp = mean_asp_rat_counts # This was 
        # probabilities_rot = mean_rot_ang_counts

        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$ \\hat{P} = %.3f, \\hat{\\gamma} = %.3f , \\hat{\\omega} = %.3f \$", pressure_value, γ_value, ω_value)
        if plot
            mat"""
            mean_asp_rat_counts = $(probabilities_asp);
            mean_asp_rat_bins = $(plot_bins_asp);
            mean_rot_ang_counts = $(probabilities_rot);
            mean_rot_ang_bins = $(plot_bins_rot);
            marker_color = $(marker_color);
            pressure_label = $(pressure_label);
            
            figure(figure_aspect_ratio);
            plot(mean_asp_rat_bins, mean_asp_rat_counts, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);

            figure(figure_rotation_angle);
            plot(mean_rot_ang_bins, mean_rot_ang_counts, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            """
        end
    end
    if plot
        mat"""
        figure(figure_aspect_ratio);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');

        figure(figure_rotation_angle);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return probabilities_asp , probabilities_rot, plot_bins_asp, plot_bins_rot
end

function plot_ellipse_ωγ_2d(γ_value) 

    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
    plot_gamma = γ_value
    gamma_value = γ_value

    # Get a list of unique input pressures
    pressure_list = unique([entry.pressure for entry in matching_γ_data]) # goes through each entry of simulation_data and get the P value at that entry
    

    # Limit range to data
    upper_limit_line_x = [1*γ_value; 1*γ_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*γ_value; .1*γ_value]
    lower_limit_line_y = [1E-5; 1]

    # Start MATLAB session
    mat"""
    figure_main = figure;
    tiled_main = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'none'); % 3 rows, 1 column

    % Axes for Attenuation
    ax_attenuation = nexttile;
    hold(ax_attenuation, 'on');
    % xlabel(ax_attenuation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel(ax_attenuation, '\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}}\$', "FontSize", 20, "Interpreter", "latex");
    set(ax_attenuation, 'XScale', 'log');
    set(ax_attenuation, 'YScale', 'log')
    set(get(ax_attenuation, 'ylabel'), 'rotation', 0);
    grid(ax_attenuation, 'on');
    box(ax_attenuation, 'on');
    %plot(ax_attenuation, $(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$');
    %plot(ax_attenuation, $(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$');
    set(ax_attenuation, 'XTickLabel', []);

    % Axes for Rotation Angle
    ax_rotation = nexttile;
    hold(ax_rotation, 'on');
    % xlabel(ax_rotation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel(ax_rotation, '\$ \\left|\\overline{\\theta} \\right|\$', "FontSize", 20, "Interpreter", "latex");
    set(ax_rotation, 'XScale', 'log');
    set(get(ax_rotation, 'ylabel'), 'rotation', 0);
    grid(ax_rotation, 'on');
    box(ax_rotation, 'on');
    set(ax_rotation, 'XTickLabel', []);

    % Axes for Aspect Ratio
    ax_aspect_ratio = nexttile;
    hold(ax_aspect_ratio, 'on');
    xlabel(ax_aspect_ratio, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel(ax_aspect_ratio, '\$ \\overline{\\frac{b}{a}} \$', "FontSize", 20, "Interpreter", "latex");
    set(ax_aspect_ratio, 'XScale', 'log');
    set(get(ax_aspect_ratio, 'ylabel'), 'rotation', 0);
    grid(ax_aspect_ratio, 'on');
    box(ax_aspect_ratio, 'on');
    """

    # get a range for plotting color from 0 to 1
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_aspect_ratio_list = Float64[];
        loop_mean_rotation_angles = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))

        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current pressure value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            jvalue_mean_aspect_ratio = mean(entry.mean_aspect_ratio for entry in matching_omega_gamma_data)
            jvalue_mean_rotation_angle = mean(entry.mean_rotation_angles for entry in matching_omega_gamma_data)
            jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)

            
            # Append values using push!
            push!(loop_mean_aspect_ratio_list, jvalue_mean_aspect_ratio)
            push!(loop_mean_rotation_angles, jvalue_mean_rotation_angle)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        loop_mean_aspect_ratio_list = loop_mean_aspect_ratio_list[valid_indices]
        loop_mean_rotation_angles = loop_mean_rotation_angles[valid_indices]
        loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]

        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$\\hat{P} = %.3f, \\hat{\\gamma} = %.2f\$", pressure_value, plot_gamma)

        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        mean_aspect_ratio = $(loop_mean_aspect_ratio_list);
        mean_rotation_angles = $(loop_mean_rotation_angles);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);

        % Plot Attenuation
        loglog(ax_attenuation, omega_gamma, mean_attenuation_x, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
        % Plot Rotation Angle
        plot(ax_rotation, omega_gamma, mean_rotation_angles, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
        % Plot Aspect Ratio
        plot(ax_aspect_ratio, omega_gamma, mean_aspect_ratio, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        """
    end

    # Add legends to the plots
    mat"""
    legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    %legend(ax_rotation, 'show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    %legend(ax_aspect_ratio, 'show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end

function plot_ωγ_attenuation_2d(gamma_value; plot=true, simulation_data=simulation_data)  # Need to fix the legend
    # Initialize outputs
    matching_omega_gamma_list = []
    loop_mean_attenuation_list = []

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
            loop_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)

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

function plot_ωγ_wavespeed_2d(gamma_value) # Need to fix lgend

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

function plot_ellipse_width_effect(ω_value, γ_value)
    width_list = [10, 20 , 50]

    normalized_variable = (width_list .- minimum(width_list)) ./ (maximum(width_list) .- minimum(width_list))

    # Intialized the plots to loop over
    mat"""
    figure_aspect_ratio = figure;
    xlabel('\$ b/a \$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ p(b/a) \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'YScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    hold on

    figure_rotation_angle = figure;
    xlabel('\$ \\theta \$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ p(\\theta) \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'YScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    hold on
    """
    
    for width in width_list

        # Assign a color
        idx = findfirst(idx -> idx ==width, width_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]
        
        simulation_data = load_data("out/processed/2d_bi_K100_W$(width).jld2")
        probabilities_asp , probabilities_rot, plot_bins_asp, plot_bins_rot = plot_ellipse_pdf(ω_value, γ_value, plot=false, simulation_data = simulation_data) 

        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\textrm{Width} = %.3f, \\hat{\\gamma} = %.3f , \\hat{\\omega} = %.3f \$", width, γ_value, ω_value)

        mat"""
        mean_asp_rat_counts = $(probabilities_asp);
        mean_asp_rat_bins = $(plot_bins_asp);
        mean_rot_ang_counts = $(probabilities_rot);
        mean_rot_ang_bins = $(plot_bins_rot);
        marker_color = $(marker_color);
        legend_label = $(legend_label);
        
        figure(figure_aspect_ratio);
        plot(mean_asp_rat_bins, mean_asp_rat_counts, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);

        figure(figure_rotation_angle);
        plot(mean_rot_ang_bins, mean_rot_ang_counts, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
        """
    end

    mat"""
    figure(figure_aspect_ratio);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');

    figure(figure_rotation_angle);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end

# Construction Zone ------------------


function plot_attenuation_width_effect(γ_value)
    width_list = [10, 20 , 50]

    normalized_variable = (width_list .- minimum(width_list)) ./ (maximum(width_list) .- minimum(width_list))

    # Intialized the plots to loop over
    mat"""
    figure_attenuation = figure;
    xlabel('\$ \\hat{\\omega}\\hat{\\gamma} \$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\frac{\\alpha}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log')
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    hold on

    """
    
    for width in width_list

        # Assign a color
        idx = findfirst(idx -> idx ==width, width_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]
        
        simulation_data = load_data("out/processed/2d_bi_K100_W$(width).jld2")
        matching_omega_gamma_list, loop_mean_attenuation_list = plot_ωγ_attenuation_2d(γ_value, plot=false, simulation_data=simulation_data)
        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\textrm{Width} = %.3f, \\hat{\\gamma} = %.3f , \\hat{\\omega} = %.3f \$", width, γ_value)

        mat"""
        matching_omega_gamma_list = $(matching_omega_gamma_list);
        loop_mean_attenuation_list = $(loop_mean_attenuation_list);
        marker_color = $(marker_color);
        legend_label = $(legend_label);
        
        figure(figure_attenuation);
        plot(matching_omega_gamma_list, loop_mean_attenuation_list, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
        """
    end

    mat"""
    figure(figure_attenuation);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end


function ML_train_ellispe()
    # CLEAN THE DATA
    # Remove rows with missing values
    dataframe_clean = filter(row -> all(x -> !(x isa Number && isnan(x)), row), data_frame) 
    dataframe_clean = dropmissing(dataframe_clean)

    # SPLIT THE DATA
    # choose the data that will be the features and targets
    feature_matrix, target_vector = unpack(dataframe_clean, ==(:mean_aspect_ratio), colname -> colname in [:mean_rotation_angles])
    feature_matrix = DataFrame(mean_rotation_angles = feature_matrix)

    # Add a new label column based on the mean_aspect_ratio
    dataframe_clean.label = ifelse.(dataframe_clean.mean_aspect_ratio .> 0.5, "high", "low")

    # NORMALIZE DATA
    standardizer = Standardizer()
    fit!(feature_matrix)
    feature_matrix = MLJModels.fit_transform!(standardizer, feature_matrix)



end


