module GranMA

using DataFrames
using CSV
using GLMakie
using LaTeXStrings
using Debugger # REPL: Debugger.@run function(); @bp
using MATLAB
using Statistics
using Printf

export plot_ellipse_ωγ_2d, random_aspect_ratio_check_2d, simulation_2d, plot_ωγ_attenuation_2d, plot_ωγ_wavespeed_2d, pack_poly_2d

"""
Creates a 2D polydisperse packing of particles with diameters ranging uniformly from `D` to `D + G`.

## Arguments
- `N::Int`: The number of particles to simulate.
- `K::Int`: The spring constant used in the simulation.
- `D::Float64`: The minimum diameter of the particles.
- `G::Float64`: The range added to `D` to determine the maximum particle diameter (`D + G`).
- `M::Float64`: The mass of each particle.
- `P_target::Float64`: The target pressure for the system.
- `W_factor::Float64`: The width factor of the simulation area.
- `seed::Int`: The random seed for reproducibility of the simulation.
- `plotit::Bool`: Flag indicating whether to plot the results (`true` for plotting, `false` otherwise).

## Description
The `pack_poly_2d` function sets up a 2D simulation environment for polydisperse particle packing. 
It uses a uniform distribution of particle diameters within the specified range and applies physical 
properties like spring constant and mass to simulate interactions. 

## Example
```julia
pack_poly_2d(5000, 100, 1.0, 0.5, 1.0, 0.01, 5.0, 1, true)
"""
function pack_poly_2d(N, K, D, G, M, P_target, W_factor, seed, plotit)
    mat"""
    addpath('src/matlab_functions/')
    N =double($(N))
    K = double($(K))
    M = double($(M))
    D = double($(D))
    G = double($(G))
    P_target = double($(P_target))
    W_factor = double($(W_factor))
    seed = double($(seed))
    plotit = double($(plotit))
    packing_poly_2d(N, K, D, G, M, P_target, W_factor, seed, plotit)
    """
end

"""
Runs a 2D simulation using MATLAB functions.

# Arguments
- `K::Int`: Spring constant.
- `M::Int`: Mass.
- `Bv::Float64`: Viscous damping coefficient.
- `w_D::Int`: Driving frequency.
- `N::Int`: Number of particles.
- `P::Float64`: Pressure parameter.
- `W::Int`: Width of the simulation area.
- `seed::Int`: Random seed for reproducibility.

# Example Function
simulation_2d(100, 1, .5, 1, 5000, .01, 5, 1)

# NOTE
Need to ensure you have the corresponding input files in the ./in/ directory 
"""
function simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    mat"""
    addpath('src/matlab_functions/')
    K = double($(K))
    M = double($(M))
    Bv = double($(Bv))
    w_D = double($(w_D))
    N =double($(N))
    P = double($(P))
    W = double($(W))
    seed = double($(seed))
    simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    """
end


"""
Plots ωγ vs Ellipse Data

# Arguments
- `data_frame::DataFrame`: This is the output from the "process_outputs_2d.m" function You will need to load the data using "data_frame = CSV.read("out/processed_simulation/your_file_name.csv", DataFrame)
.
- `gamma_value::Float65`: γ value that you want to plot.
- `flag::Any`: Add a third flag if you want to use MATLAB to plot.

# Example Function
plot_ellipse_ωγ_2d(data_frame, .1)

# NOTE
Need to ensure you have the corresponding input files in the ./in/ directory 
"""
function plot_ellipse_ωγ_2d(data_frame, gamma_value, flag::Any) # using MATLAB

    # Define parameters to plot
    pressure_list = sort(unique(data_frame.input_pressure))
    plot_pressure = pressure_list
    gamma_list = sort(unique(data_frame.gamma))
    closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
    plot_gamma = gamma_list[closest_gamma_match_index]

    # Filter the table to only those data
    matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
    combined_index = matching_gamma_index
    filtered_data_frame = data_frame[combined_index, :]

    # Limit range to data
    upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
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

    # Normalize the gamma values
    normalized_variable = (log.(plot_pressure) .- minimum(log.(plot_pressure))) ./ (maximum(log.(plot_pressure)) .- minimum(log.(plot_pressure)))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # For idx, only show current pressure data
        iloop_pressure_value = plot_pressure[idx]
        iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
        iloop_combined_index = iloop_pressure_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_aspect_ratio_list = Float64[];
        loop_mean_rotation_angles = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

        for jdx in eachindex(iloop_omega_gamma_list)

            # Get the idex for the current omega gamma value
            matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
            jloop_data_frame = iloop_data_frame[matching_jdx,:]

            # get the mean over all seeds
            jvalue_mean_aspect_ratio = mean(jloop_data_frame.mean_aspect_ratio)
            jvalue_mean_rotation_angle = mean(jloop_data_frame.mean_rotation_angles)
            jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)
            
            # Append values using push!
            push!(loop_mean_aspect_ratio_list, jvalue_mean_aspect_ratio)
            push!(loop_mean_rotation_angles, jvalue_mean_rotation_angle)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        valid_indices = iloop_omega_gamma_list .<= gamma_value.*2
        iloop_omega_gamma_list = iloop_omega_gamma_list[valid_indices]
        loop_mean_aspect_ratio_list = loop_mean_aspect_ratio_list[valid_indices]
        loop_mean_rotation_angles = loop_mean_rotation_angles[valid_indices]
        loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]
        
        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$\\hat{P} = %.3f, \\hat{\\gamma} = %.2f\$", iloop_pressure_value, plot_gamma)

        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(iloop_omega_gamma_list);
        mean_aspect_ratio = $(loop_mean_aspect_ratio_list);
        mean_rotation_angles = $(loop_mean_rotation_angles);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(iloop_pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);

        % Plot Attenuation
        plot(ax_attenuation, omega_gamma, mean_attenuation_x, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
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


"""
Picks 5 random lines within a pre-loaded data frame that only vary by seed.

# Arguments
- None.

# Example Function
random_aspect_ratio_check_2d()

# NOTE
Need to ensure you have the corresponding input files in the ./in/ directory 
"""
function random_aspect_ratio_check_2d()
    # pick a random row 
    random_row = data_frame[rand(1:nrow(data_frame)) ,:]

    # seeds have input_pressure and omega gamma in common
    common_rows = data_frame[in.(data_frame.input_pressure, random_row.input_pressure) .& in.(data_frame.omegagamma, random_row.omegagamma),:]

    return common_rows    
end

"""
Plots ωγ vs α / ω

# Arguments
- `data_frame::DataFrame`: This is the output from the "process_outputs_2d.m" function You will need to load the data using "data_frame = CSV.read("out/processed_simulation/your_file_name.csv", DataFrame)
- `gamma_value::Float65`: γ value that you want to plot.
- `flag::Any`: Add a third flag if you want to use MATLAB to plot.

# Example Function
plot_ωγ_attenuation_2d(data_frame, .1)

# NOTE
Need to ensure you have the corresponding input files in the ./in/ directory 
"""
function plot_ωγ_attenuation_2d(data_frame, gamma_value)

    # Define parameters to plot
    pressure_list = sort(unique(data_frame.input_pressure))
    plot_pressure = pressure_list
    gamma_list = sort(unique(data_frame.gamma))
    closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
    plot_gamma = gamma_list[closest_gamma_match_index]

    # Filter the table to only those data
    matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
    combined_index = matching_gamma_index
    filtered_data_frame = data_frame[combined_index, :]

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);

    fig = Figure()
    ax = Axis(fig[1,1], 
        xlabel=L"\hat{\omega}\hat{\gamma}",
        ylabel=L"\frac{\hat{\alpha}}{\hat{\omega}}", 
        xscale = log10,
        yscale = log10,
        ylabelrotation=0,
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(10),
        yminorticksvisible = true,
        yminorgridvisible = true,
        yminorticks = IntervalsBetween(10))
    lines!(ax,theory_x, theory_y)

    # Normalize the gamma values
    normalized_variable = (plot_pressure .- minimum(plot_pressure)) ./ (maximum(plot_pressure) .- minimum(plot_pressure))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = RGBf(normalized_variable[idx], 0, 1-normalized_variable[idx]);

        # For idx, only show current pressure data
        iloop_pressure_value = plot_pressure[idx]
        iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
        iloop_combined_index = iloop_pressure_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];
        loop_omegagamma_plot_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

        for jdx in eachindex(iloop_omega_gamma_list)

            # Get the idex for the current omega gamma value
            matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
            jloop_data_frame = iloop_data_frame[matching_jdx,:]

            # get the mean over all seeds
            jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)

            # Append values only if they are positive
            if jvalue_mean_alphaoveromega >= 0
                push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
                push!(loop_omegagamma_plot_list, iloop_omega_gamma_list[jdx])
            end
        end
        scatterlines!(ax, loop_omegagamma_plot_list, loop_mean_attenuation_list, color=marker_color, label=L"\hat{P} = %$(iloop_pressure_value)")
    end
    axislegend(position = :rb)
    fig

end

function plot_ωγ_attenuation_2d(data_frame, gamma_value, flag::Any) # using MATLAB

    # Define parameters to plot
    pressure_list = sort(unique(data_frame.input_pressure))
    plot_pressure = pressure_list
    gamma_list = sort(unique(data_frame.gamma))
    closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
    plot_gamma = gamma_list[closest_gamma_match_index]

    # Filter the table to only those data
    matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
    combined_index = matching_gamma_index
    filtered_data_frame = data_frame[combined_index, :]

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    lower_limit_line_y = [1E-5; 1]

    # Intialized the plots to loop over
    mat"""
    figure_attenuation = figure;
    loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory'), hold on
    plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    """

    # Normalize the gamma values
    normalized_variable = (log.(plot_pressure) .- minimum(log.(plot_pressure))) ./ (maximum(log.(plot_pressure)) .- minimum(log.(plot_pressure)))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]];

        # For idx, only show current pressure data
        iloop_pressure_value = plot_pressure[idx]
        iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
        iloop_combined_index = iloop_pressure_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

        for jdx in eachindex(iloop_omega_gamma_list)

            # Get the idex for the current omega gamma value
            matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
            jloop_data_frame = iloop_data_frame[matching_jdx,:]

            # get the mean over all seeds
            jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)
            @bp
            # Append values using push!
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end
        
        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(iloop_omega_gamma_list);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(iloop_pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color= $(marker_color);
        pressure_label = sprintf('Pressure = %.2f, Gamma = %.2f (Aspect Ratio)', $(iloop_pressure_value), $(plot_gamma));
        % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
        pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Attenuation)', $(iloop_pressure_value), $(plot_gamma));
        
        figure(figure_attenuation);
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_attenuation_x, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
        """
    end

    # Add legends to the plots
    mat"""
    % Add legends to the MATLAB plots
    figure(figure_attenuation);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end

"""
Plots ωγ vs c

# Arguments
- `data_frame::DataFrame`: This is the output from the "process_outputs_2d.m" function You will need to load the data using "data_frame = CSV.read("out/processed_simulation/your_file_name.csv", DataFrame)
.
- `gamma_value::Float65`: γ value that you want to plot.
- `flag::Any`: Add a third flag if you want to use MATLAB to plot.

# Example Function
plot_ωγ_wavespeed_2d(data_frame, .1)

# NOTE
Need to ensure you have the corresponding input files in the ./in/ directory 
"""
function plot_ωγ_wavespeed_2d(data_frame, gamma_value)

    # Define parameters to plot
    pressure_list = sort(unique(data_frame.input_pressure))
    plot_pressure = pressure_list
    gamma_list = sort(unique(data_frame.gamma))
    closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
    plot_gamma = gamma_list[closest_gamma_match_index]

    # Filter the table to only those data
    matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
    combined_index = matching_gamma_index
    filtered_data_frame = data_frame[combined_index, :]

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = 1 ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(0.5) ./ (1 .+ theory_x.^2);

    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel=L"\hat{\omega}\hat{\gamma}",
        ylabel=L"\hat{c}", 
        xscale = log10,
        yscale = log10,
        yminorticksvisible = true,
        yminorgridvisible = true,
        yminorticks = IntervalsBetween(5),
        xminorticksvisible = true,
        xminorgridvisible = true,
        xminorticks = IntervalsBetween(5))

        lines!(ax,theory_x, 1 ./theory_y)

    # Normalize the gamma values
    normalized_variable = (log.(plot_pressure) .- log.(minimum(plot_pressure))) ./ (maximum(log.(plot_pressure)) .- minimum(log.(plot_pressure)))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = RGBf(normalized_variable[idx], 0, 1-normalized_variable[idx]);

        # For idx, only show current pressure data
        iloop_pressure_value = plot_pressure[idx]
        iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
        iloop_combined_index = iloop_pressure_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_wavespeed_list = Float64[];
        loop_omegagamma_plot_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

        for jdx in eachindex(iloop_omega_gamma_list)

            # Get the idex for the current omega gamma value
            matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
            jloop_data_frame = iloop_data_frame[matching_jdx,:]

            # get the mean over all seeds
            jvalue_mean_wavespeed = mean(jloop_data_frame.wavespeed_x)

            # Append values only if they are positive
            if jvalue_mean_wavespeed >= 0
                push!(loop_mean_wavespeed_list, jvalue_mean_wavespeed)
                push!(loop_omegagamma_plot_list, iloop_omega_gamma_list[jdx])
            end
        end
        scatterlines!(ax, loop_omegagamma_plot_list, loop_mean_wavespeed_list, color=marker_color, label=L"\hat{P} = %$(iloop_pressure_value), \hat{\gamma}=%$(gamma_value)")

    end
    axislegend(position = :lt)
    fig
end

function plot_ωγ_wavespeed_2d(data_frame, gamma_value, flag::Any) # using MATLAB

    # Define parameters to plot
    pressure_list = sort(unique(data_frame.input_pressure))
    plot_pressure = pressure_list
    gamma_list = sort(unique(data_frame.gamma))
    closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
    plot_gamma = gamma_list[closest_gamma_match_index]

    # Filter the table to only those data
    matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
    combined_index = matching_gamma_index
    filtered_data_frame = data_frame[combined_index, :]

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = 1 ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(0.5) ./ (1 .+ theory_x.^2);
    upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    upper_limit_line_y = [9E-1; 2]
    lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    lower_limit_line_y = [9E-1; 2]

    # Intialized the plots to loop over
    mat"""
    figure_wavespeed = figure;
    loglog($(theory_x), 1./$(theory_y), 'k', 'DisplayName', '1-D Theory'), hold on
    plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$\\hat{c} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    """

    # Normalize the gamma values
    normalized_variable = (plot_pressure .- minimum(plot_pressure)) ./ (maximum(plot_pressure) .- minimum(plot_pressure))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]];

        # For idx, only show current pressure data
        iloop_pressure_value = plot_pressure[idx]
        iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
        iloop_combined_index = iloop_pressure_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_wavespeed_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

        for jdx in eachindex(iloop_omega_gamma_list)

            # Get the idex for the current omega gamma value
            matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
            jloop_data_frame = iloop_data_frame[matching_jdx,:]

            # get the mean over all seeds
            jvalue_mean_wavespeed = mean(jloop_data_frame.wavespeed_x)
            @bp
            # Append values using push!
            push!(loop_mean_wavespeed_list, jvalue_mean_wavespeed)
        end
        
        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(iloop_omega_gamma_list);
        mean_wavespeed_x = $(loop_mean_wavespeed_list);
        iloop_pressure_value = $(iloop_pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = sprintf('Wavespeed X = %.2f, Gamma = %.2f (Aspect Ratio)', $(iloop_pressure_value), $(plot_gamma));
        % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
        pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Wavespeed X)', $(iloop_pressure_value), $(plot_gamma));
        
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

function plot_collapsed_ellipse_ωγ_2d(data_frame, gamma_value, flag::Any) # using MATLAB

    # Define parameters to plot
    pressure_list = sort(unique(data_frame.input_pressure))
    plot_pressure = pressure_list
    gamma_list = sort(unique(data_frame.gamma))
    closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
    plot_gamma = gamma_list[closest_gamma_match_index]

    # Filter the table to only those data
    matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
    combined_index = matching_gamma_index
    filtered_data_frame = data_frame[combined_index, :]

    # Limit range to data
    upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    lower_limit_line_y = [1E-5; 1]

    # Start MATLAB session
    mat"""
    % Initialize the MATLAB figure for aspect ratio
    figure_aspect_ratio = figure;
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\overline{\\frac{b}{a}} \$', "FontSize", 20, "Interpreter", "latex");
    plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;

    % Initialize the MATLAB figure for rotation angle
    figure_rotation_angle = figure;
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\left|\\overline{\\theta} \\right|\$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    """

    # Normalize the gamma values
    normalized_variable = (log.(plot_pressure) .- minimum(log.(plot_pressure))) ./ (maximum(log.(plot_pressure)) .- minimum(log.(plot_pressure)))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # For idx, only show current pressure data
        iloop_pressure_value = plot_pressure[idx]
        iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
        iloop_combined_index = iloop_pressure_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_aspect_ratio_list = Float64[];
        loop_mean_rotation_angles = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

        for jdx in eachindex(iloop_omega_gamma_list)

            # Get the idex for the current omega gamma value
            matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
            jloop_data_frame = iloop_data_frame[matching_jdx,:]

            # get the mean over all seeds
            jvalue_mean_aspect_ratio = mean(jloop_data_frame.mean_aspect_ratio)
            jvalue_mean_rotation_angle = mean(jloop_data_frame.mean_rotation_angles)
            jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)
            @bp
            # Append values using push!
            push!(loop_mean_aspect_ratio_list, jvalue_mean_aspect_ratio)
            push!(loop_mean_rotation_angles, jvalue_mean_rotation_angle)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end
        
        # Transfer data to MATLAB
        mat"""
        % Transfer data to MATLAB
        omega_gamma = $(iloop_omega_gamma_list);
        mean_aspect_ratio = $(loop_mean_aspect_ratio_list);
        mean_rotation_angles = $(loop_mean_rotation_angles);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(iloop_pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = sprintf('Pressure = %.2f, Gamma = %.2f (Aspect Ratio)', $(iloop_pressure_value), $(plot_gamma));
        % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
        pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Attenuation)', $(iloop_pressure_value), $(plot_gamma));

        
        figure(figure_aspect_ratio);
        yyaxis left;
        plot(omega_gamma, mean_aspect_ratio, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        yyaxis right;
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_attenuation_x, '--x','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
        
        figure(figure_rotation_angle);
        yyaxis left;
        plot(omega_gamma, mean_rotation_angles, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        yyaxis right;
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_attenuation_x, '--x','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
        """
    end

    # Add legends to the plots
    mat"""
    % Add legends to the MATLAB plots
    figure(figure_aspect_ratio);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');

    figure(figure_rotation_angle);
    plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end
end