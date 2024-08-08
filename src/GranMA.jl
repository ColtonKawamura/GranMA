module GranMA

using DataFrames
using CSV
# using GLMakie # Not supported on HPC. Need to comment out before running on HPC
using LaTeXStrings
using Debugger # REPL: Debugger.@run function(); @bp
using MATLAB
using Statistics
using Printf
using MAT
using Glob
using JLD2
using IterTools

export  generate_simulation_jobs,
        load_data2,
        load_data,
        file_data,
        crunch_and_save,
        crunch,
        save_data,
        simulation_2d,
        pack_poly_2d, 
        plot_ellipse_low_pressure,
        process_outputs_2d

mutable struct file_data
    pressure::Float64
    omega::Float64
    gamma::Float64
    asp_rat_counts::Vector{Float64}
    asp_rat_bins::Vector{Float64}
    rot_ang_counts::Vector{Float64}
    rot_ang_bins::Vector{Float64}
    omega_gamma::Float64
    seed::Float64
    pressure_actual::Float64
    attenuation_x::Float64
    attenuation_y::Float64
    wavespeed_x::Float64
    wavenumber_x::Float64
    mean_aspect_ratio::Float64
    mean_rotation_angles::Float64
    fft_limit_x::Float64
    ellipse_stats::Matrix{Float64}
    fft_limit_y::Float64
    wavenumber_y::Float64
    alphaoveromega_x::Float64
    alphaoveromega_y::Float64
end

function generate_simulation_jobs(filename::String, K_values::Vector{T1}, M_values::Vector{T2}, Bv_values::Vector{T3}, w_D_values::Vector{T4}, N_values::Vector{T5}, P_values::Vector{T6}, W_values::Vector{T7}, seeds::Vector{T8}) where {T1, T2, T3, T4, T5, T6, T7, T8}
    # generate_simulation_jobs("testfile.txt", [100],[1],exp10.(-4:.5:1),exp10.(-4:.5:1),[10000],[.1],[10,20,50,100],[1,2,3,4,5])
    # Convert input vectors to the required types
    K_values = Int64.(K_values)
    M_values = Int64.(M_values)
    Bv_values = Float64.(Bv_values)
    w_D_values = Float64.(w_D_values)
    N_values = Float64.(N_values)
    P_values = Float64.(P_values)
    W_values = Int64.(W_values)
    seeds = Int64.(seeds)

    # Function to generate the MATLAB command
    function generate_matlab_command(K, M, Bv, w_D, N, P, W, seed)
        return "matlab -nodisplay -nosplash -r \"addpath('./src/matlab_functions/'); packing_bi_2d($K, $M, $Bv, $w_D, $N, $P, $W, $seed, false); exit\""
    end

    # Generate all combinations of parameters using IterTools.product
    combinations = IterTools.product(K_values, M_values, Bv_values, w_D_values, N_values, P_values, W_values, seeds)

    # Open a file to write the MATLAB commands
    open(filename, "w") do file
        for combo in combinations
            command = generate_matlab_command(combo...)
            println(file, command)
        end
    end

    println("Commands written to $filename")
end

function crunch_and_save()
    K = 100
    simulation_data = crunch()
    save_data(simulation_data, K)
    
    # Load the data back to verify
    reloaded_data = load_data(K)
    println("Data reloaded successfully. Number of entries: ", length(reloaded_data))
end

function crunch()
    directory = "out/simulation_2d/bi_K100_all_seeds/"
    # directory = "out/simulation_2d/"
    mat_files = glob("*.mat", directory)
    simulation_data = file_data[]
    
    for file_name in mat_files
        iloop_file_data = matread(file_name)

        # Extract data fields
        input_pressure = iloop_file_data["input_pressure"]
        omega = iloop_file_data["driving_angular_frequency_dimensionless"]
        gamma = iloop_file_data["gamma_dimensionless"]
        pressure_actual = iloop_file_data["pressure_dimensionless"]
        attenuation_x = iloop_file_data["attenuation_x_dimensionless"]
        attenuation_y = iloop_file_data["attenuation_y_dimensionless"]
        wavespeed_x = iloop_file_data["wavespeed_x"]
        wavenumber_x = iloop_file_data["wavenumber_x_dimensionless"]
        mean_aspect_ratio = iloop_file_data["mean_aspect_ratio"]
        mean_rotation_angles = iloop_file_data["mean_rotation_angles"]
        wavenumber_y = iloop_file_data["wavenumber_y_dimensionless"]

        # Check for NaN values
        if isnan(mean_rotation_angles) || isnan(mean_aspect_ratio) || 
           isnan(wavenumber_x) || isnan(wavenumber_y) ||
           isnan(wavespeed_x) || isnan(input_pressure) || 
           isnan(omega) || isnan(gamma) || isnan(pressure_actual) || 
           isnan(attenuation_x) || isnan(attenuation_y)
            println("Skipping file due to NaN values: $file_name")
            continue
        end

        # Convert matrix-type data to vectors
        asp_rat_counts = vec(iloop_file_data["asp_rat_counts"])
        asp_rat_bins = vec(iloop_file_data["asp_rat_bins"])
        rot_ang_counts = vec(iloop_file_data["rot_ang_counts"])
        rot_ang_bins = vec(iloop_file_data["rot_ang_bins"])

        # Handle potentially empty arrays for fft limits
        fft_x = get(iloop_file_data, "initial_distance_from_oscillation_output_x_fft", [])
        fft_y = get(iloop_file_data, "initial_distance_from_oscillation_output_y_fft", [])
        fft_limit_x = isempty(fft_x) ? NaN : maximum(fft_x)
        fft_limit_y = isempty(fft_y) ? NaN : maximum(fft_y)

        data_entry = file_data(
            input_pressure,
            omega,
            gamma,
            asp_rat_counts,
            asp_rat_bins,
            rot_ang_counts,
            rot_ang_bins,
            omega * gamma,
            iloop_file_data["seed"],
            pressure_actual,
            -attenuation_x,
            -attenuation_y,
            -wavespeed_x,
            -wavenumber_x,
            mean_aspect_ratio,
            mean_rotation_angles,
            fft_limit_x,
            iloop_file_data["ellipse_stats_nonzero"],
            fft_limit_y,
            wavenumber_y,
            -attenuation_x / omega,
            -attenuation_y / omega
        )

        push!(simulation_data, data_entry)
    end

    return simulation_data
end

function save_data(simulation_data::Vector{file_data}, K::Int)
    file_name = "out/processed/2d_K$(K).jld2"
    @save file_name simulation_data
    println("Data saved to $file_name")
end

function load_data(K::Int)::Vector{file_data}
    file_name = "out/processed/2d_K$(K).jld2"
    @load file_name simulation_data
    return simulation_data
end

# function load_data2(K::Int)::Vector{file_data}
#     file_name = "out/processed/2d_K$(K).jld2"
#     simulation_data = JLD2.load(file_name, "simulation_data")
    
#     # Convert to the correct type
#     return [file_data(d.pressure, d.omega, d.gamma, d.asp_rat_counts, d.asp_rat_bins,
#                       d.rot_ang_counts, d.rot_ang_bins, d.omega_gamma, d.seed,
#                       d.pressure_actual, d.attenuation_x, d.attenuation_y, d.wavespeed_x,
#                       d.wavenumber_x, d.mean_aspect_ratio, d.mean_rotation_angles,
#                       d.fft_limit_x, d.ellipse_stats, d.fft_limit_y, d.wavenumber_y)
#             for d in simulation_data]
# end



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


# function plot_ellipse_ωγ_2d(data_frame, gamma_value) # using MATLAB

#     # Define parameters to plot
#     pressure_list = sort(unique(data_frame.input_pressure))
#     plot_pressure = pressure_list
#     gamma_list = sort(unique(data_frame.gamma))
#     closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
#     plot_gamma = gamma_list[closest_gamma_match_index]

#     # Filter the table to only those data
#     matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
#     combined_index = matching_gamma_index
#     filtered_data_frame = data_frame[combined_index, :]

#     # Limit range to data
#     upper_limit_line_x = [1*gamma_value; 1*gamma_value]
#     upper_limit_line_y = [1E-5; 1]
#     lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
#     lower_limit_line_y = [1E-5; 1]

#     # Start MATLAB session
#     mat"""
#     figure_main = figure;
#     tiled_main = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'none'); % 3 rows, 1 column

#     % Axes for Attenuation
#     ax_attenuation = nexttile;
#     hold(ax_attenuation, 'on');
#     % xlabel(ax_attenuation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
#     ylabel(ax_attenuation, '\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}}\$', "FontSize", 20, "Interpreter", "latex");
#     set(ax_attenuation, 'XScale', 'log');
#     set(ax_attenuation, 'YScale', 'log')
#     set(get(ax_attenuation, 'ylabel'), 'rotation', 0);
#     grid(ax_attenuation, 'on');
#     box(ax_attenuation, 'on');
#     %plot(ax_attenuation, $(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$');
#     %plot(ax_attenuation, $(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$');
#     set(ax_attenuation, 'XTickLabel', []);

#     % Axes for Rotation Angle
#     ax_rotation = nexttile;
#     hold(ax_rotation, 'on');
#     % xlabel(ax_rotation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
#     ylabel(ax_rotation, '\$ \\left|\\overline{\\theta} \\right|\$', "FontSize", 20, "Interpreter", "latex");
#     set(ax_rotation, 'XScale', 'log');
#     set(get(ax_rotation, 'ylabel'), 'rotation', 0);
#     grid(ax_rotation, 'on');
#     box(ax_rotation, 'on');
#     set(ax_rotation, 'XTickLabel', []);

#     % Axes for Aspect Ratio
#     ax_aspect_ratio = nexttile;
#     hold(ax_aspect_ratio, 'on');
#     xlabel(ax_aspect_ratio, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
#     ylabel(ax_aspect_ratio, '\$ \\overline{\\frac{b}{a}} \$', "FontSize", 20, "Interpreter", "latex");
#     set(ax_aspect_ratio, 'XScale', 'log');
#     set(get(ax_aspect_ratio, 'ylabel'), 'rotation', 0);
#     grid(ax_aspect_ratio, 'on');
#     box(ax_aspect_ratio, 'on');
#     """

#     # Normalize the gamma values
#     normalized_variable = (log.(plot_pressure) .- minimum(log.(plot_pressure))) ./ (maximum(log.(plot_pressure)) .- minimum(log.(plot_pressure)))

#     # Create a line for each gamma value across all pressure_list
#     for idx in eachindex(plot_pressure)

#         marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

#         # For idx, only show current pressure data
#         iloop_pressure_value = plot_pressure[idx]
#         iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
#         iloop_combined_index = iloop_pressure_index
#         iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

#         # Initizalized vectors for just this pressure
#         loop_mean_aspect_ratio_list = Float64[];
#         loop_mean_rotation_angles = Float64[];
#         loop_mean_attenuation_list = Float64[];

#         # Look at a single omega gamma value since each one spans all seeds
#         iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

#         for jdx in eachindex(iloop_omega_gamma_list)

#             # Get the idex for the current omega gamma value
#             matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
#             jloop_data_frame = iloop_data_frame[matching_jdx,:]

#             # get the mean over all seeds
#             jvalue_mean_aspect_ratio = mean(jloop_data_frame.mean_aspect_ratio)
#             jvalue_mean_rotation_angle = mean(jloop_data_frame.mean_rotation_angles)
#             jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)
            
#             # Append values using push!
#             push!(loop_mean_aspect_ratio_list, jvalue_mean_aspect_ratio)
#             push!(loop_mean_rotation_angles, jvalue_mean_rotation_angle)
#             push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
#         end

#         # Filter data to include only points where omega_gamma <= gamma_value
#         valid_indices = iloop_omega_gamma_list .<= gamma_value.*2
#         iloop_omega_gamma_list = iloop_omega_gamma_list[valid_indices]
#         loop_mean_aspect_ratio_list = loop_mean_aspect_ratio_list[valid_indices]
#         loop_mean_rotation_angles = loop_mean_rotation_angles[valid_indices]
#         loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]
        
#         # This is needed because MATLAB.jl has a hard time escaping \'s
#         pressure_label = @sprintf("\$\\hat{P} = %.3f, \\hat{\\gamma} = %.2f\$", iloop_pressure_value, plot_gamma)

#         # Transfer data to MATLAB
#         mat"""
#         omega_gamma = $(iloop_omega_gamma_list);
#         mean_aspect_ratio = $(loop_mean_aspect_ratio_list);
#         mean_rotation_angles = $(loop_mean_rotation_angles);
#         mean_attenuation_x = $(loop_mean_attenuation_list);
#         iloop_pressure_value = $(iloop_pressure_value);
#         plot_gamma = $(plot_gamma);
#         marker_color = $(marker_color);
#         pressure_label = $(pressure_label);

#         % Plot Attenuation
#         loglog(ax_attenuation, omega_gamma, mean_attenuation_x, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
#         % Plot Rotation Angle
#         plot(ax_rotation, omega_gamma, mean_rotation_angles, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
#         % Plot Aspect Ratio
#         plot(ax_aspect_ratio, omega_gamma, mean_aspect_ratio, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
#         """
#     end

#     # Add legends to the plots
#     mat"""
#     legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
#     %legend(ax_rotation, 'show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
#     %legend(ax_aspect_ratio, 'show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
#     """
# end

# function plot_ellipse_ωγ_2d(data_frame, gamma_value, flag::Any) # using MATLAB

#     # Define parameters to plot
#     pressure_list = sort(unique(data_frame.input_pressure))
#     plot_pressure = pressure_list
#     gamma_list = sort(unique(data_frame.gamma))
#     closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
#     plot_gamma = gamma_list[closest_gamma_match_index]

#     # Filter the table to only those data
#     matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
#     combined_index = matching_gamma_index
#     filtered_data_frame = data_frame[combined_index, :]

#     # Limit range to data
#     upper_limit_line_x = [1*gamma_value; 1*gamma_value]
#     upper_limit_line_y = [1E-5; 1]
#     lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
#     lower_limit_line_y = [1E-5; 1]

#     # Start MATLAB session
#     mat"""
#     % Initialize the MATLAB figure for aspect ratio
#     figure_aspect_ratio = figure;
#     hold on;
#     xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
#     ylabel('\$ \\overline{\\frac{b}{a}} \$', "FontSize", 20, "Interpreter", "latex");
#     plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
#     plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')
#     set(gca, 'XScale', 'log');
#     set(get(gca, 'ylabel'), 'rotation', 0);
#     grid on;
#     box on;

#     % Initialize the MATLAB figure for rotation angle
#     figure_rotation_angle = figure;
#     hold on;
#     xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
#     ylabel('\$\\overline{  \\left| \\theta \\right|}\$', "FontSize", 20, "Interpreter", "latex");
#     set(gca, 'XScale', 'log');
#     set(get(gca, 'ylabel'), 'rotation', 0);
#     grid on;
#     box on;
#     """

#     # Normalize the gamma values
#     normalized_variable = (log.(plot_pressure) .- minimum(log.(plot_pressure))) ./ (maximum(log.(plot_pressure)) .- minimum(log.(plot_pressure)))

#     # Create a line for each gamma value across all pressure_list
#     for idx in eachindex(plot_pressure)

#         marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

#         # For idx, only show current pressure data
#         iloop_pressure_value = plot_pressure[idx]
#         iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
#         iloop_combined_index = iloop_pressure_index
#         iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

#         # Initizalized vectors for just this pressure
#         loop_mean_aspect_ratio_list = Float64[];
#         loop_mean_rotation_angles = Float64[];
#         loop_mean_attenuation_list = Float64[];

#         # Look at a single omega gamma value since each one spans all seeds
#         iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

#         for jdx in eachindex(iloop_omega_gamma_list)

#             # Get the idex for the current omega gamma value
#             matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
#             jloop_data_frame = iloop_data_frame[matching_jdx,:]

#             # get the mean over all seeds
#             jvalue_mean_aspect_ratio = mean(jloop_data_frame.mean_aspect_ratio)
#             jvalue_mean_rotation_angle = mean(jloop_data_frame.mean_rotation_angles)
#             jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)
#             @bp
#             # Append values using push!
#             push!(loop_mean_aspect_ratio_list, jvalue_mean_aspect_ratio)
#             push!(loop_mean_rotation_angles, jvalue_mean_rotation_angle)
#             push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
#         end
        
#         # Transfer data to MATLAB
#         mat"""
#         % Transfer data to MATLAB
#         omega_gamma = $(iloop_omega_gamma_list);
#         mean_aspect_ratio = $(loop_mean_aspect_ratio_list);
#         mean_rotation_angles = $(loop_mean_rotation_angles);
#         mean_attenuation_x = $(loop_mean_attenuation_list);
#         iloop_pressure_value = $(iloop_pressure_value);
#         plot_gamma = $(plot_gamma);
#         marker_color = $(marker_color);
#         pressure_label = sprintf('Pressure = %.2f, Gamma = %.2f (Aspect Ratio)', $(iloop_pressure_value), $(plot_gamma));
#         % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
#         pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Attenuation)', $(iloop_pressure_value), $(plot_gamma));

        
#         figure(figure_aspect_ratio);
#         yyaxis left;
#         plot(omega_gamma, mean_aspect_ratio, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
#         yyaxis right;
#         set(gca, 'Yscale', 'log');
#         plot(omega_gamma, mean_attenuation_x, '--x','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
        
#         figure(figure_rotation_angle);
#         yyaxis left;
#         plot(omega_gamma, mean_rotation_angles, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
#         yyaxis right;
#         set(gca, 'Yscale', 'log');
#         plot(omega_gamma, mean_attenuation_x, '--x','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
#         """
#     end

#     # Add legends to the plots
#     mat"""
#     % Add legends to the MATLAB plots
#     figure(figure_aspect_ratio);
#     legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');

#     figure(figure_rotation_angle);
#     plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
#     plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')
#     legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
#     """
# end


# function random_aspect_ratio_check_2d()
#     # pick a random row 
#     random_row = data_frame[rand(1:nrow(data_frame)) ,:]

#     # seeds have input_pressure and omega gamma in common
#     common_rows = data_frame[in.(data_frame.input_pressure, random_row.input_pressure) .& in.(data_frame.omegagamma, random_row.omegagamma),:]

#     return common_rows    
# end


# function plot_ωγ_attenuation_2d(data_frame, gamma_value) # This is with GLMAKIE plotting package

#     # Define parameters to plot
#     pressure_list = sort(unique(data_frame.input_pressure))
#     plot_pressure = pressure_list
#     gamma_list = sort(unique(data_frame.gamma))
#     closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
#     plot_gamma = gamma_list[closest_gamma_match_index]

#     # Filter the table to only those data
#     matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
#     combined_index = matching_gamma_index
#     filtered_data_frame = data_frame[combined_index, :]

#     # Define the plot limits to match the 1D theory plot curves
#     theory_x = collect(3E-4:1E-5:3)
#     theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);

#     fig = Figure()
#     ax = Axis(fig[1,1], 
#         xlabel=L"\hat{\omega}\hat{\gamma}",
#         ylabel=L"\frac{\hat{\alpha}}{\hat{\omega}}", 
#         xscale = log10,
#         yscale = log10,
#         ylabelrotation=0,
#         xminorticksvisible = true,
#         xminorgridvisible = true,
#         xminorticks = IntervalsBetween(10),
#         yminorticksvisible = true,
#         yminorgridvisible = true,
#         yminorticks = IntervalsBetween(10))
#     lines!(ax,theory_x, theory_y)

#     # Normalize the gamma values
#     normalized_variable = (plot_pressure .- minimum(plot_pressure)) ./ (maximum(plot_pressure) .- minimum(plot_pressure))

#     # Create a line for each gamma value across all pressure_list
#     for idx in eachindex(plot_pressure)

#         marker_color = RGBf(normalized_variable[idx], 0, 1-normalized_variable[idx]);

#         # For idx, only show current pressure data
#         iloop_pressure_value = plot_pressure[idx]
#         iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
#         iloop_combined_index = iloop_pressure_index
#         iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

#         # Initizalized vectors for just this pressure
#         loop_mean_attenuation_list = Float64[];
#         loop_omegagamma_plot_list = Float64[];

#         # Look at a single omega gamma value since each one spans all seeds
#         iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

#         for jdx in eachindex(iloop_omega_gamma_list)

#             # Get the idex for the current omega gamma value
#             matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
#             jloop_data_frame = iloop_data_frame[matching_jdx,:]

#             # get the mean over all seeds
#             jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)

#             # Append values only if they are positive
#             if jvalue_mean_alphaoveromega >= 0
#                 push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
#                 push!(loop_omegagamma_plot_list, iloop_omega_gamma_list[jdx])
#             end
#         end
#         scatterlines!(ax, loop_omegagamma_plot_list, loop_mean_attenuation_list, color=marker_color, label=L"\hat{P} = %$(iloop_pressure_value)")
#     end
#     axislegend(position = :rb)
#     fig

# end

# function plot_ωγ_wavespeed_2d(data_frame, gamma_value)

#     # Define parameters to plot
#     pressure_list = sort(unique(data_frame.input_pressure))
#     plot_pressure = pressure_list
#     gamma_list = sort(unique(data_frame.gamma))
#     closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
#     plot_gamma = gamma_list[closest_gamma_match_index]

#     # Filter the table to only those data
#     matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
#     combined_index = matching_gamma_index
#     filtered_data_frame = data_frame[combined_index, :]

#     # Define the plot limits to match the 1D theory plot curves
#     theory_x = collect(3E-4:1E-5:3)
#     theory_y = 1 ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(0.5) ./ (1 .+ theory_x.^2);

#     fig = Figure()
#     ax = Axis(fig[1,1],
#         xlabel=L"\hat{\omega}\hat{\gamma}",
#         ylabel=L"\hat{c}", 
#         xscale = log10,
#         yscale = log10,
#         yminorticksvisible = true,
#         yminorgridvisible = true,
#         yminorticks = IntervalsBetween(5),
#         xminorticksvisible = true,
#         xminorgridvisible = true,
#         xminorticks = IntervalsBetween(5))

#         lines!(ax,theory_x, 1 ./theory_y)

#     # Normalize the gamma values
#     normalized_variable = (log.(plot_pressure) .- log.(minimum(plot_pressure))) ./ (maximum(log.(plot_pressure)) .- minimum(log.(plot_pressure)))

#     # Create a line for each gamma value across all pressure_list
#     for idx in eachindex(plot_pressure)

#         marker_color = RGBf(normalized_variable[idx], 0, 1-normalized_variable[idx]);

#         # For idx, only show current pressure data
#         iloop_pressure_value = plot_pressure[idx]
#         iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
#         iloop_combined_index = iloop_pressure_index
#         iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

#         # Initizalized vectors for just this pressure
#         loop_mean_wavespeed_list = Float64[];
#         loop_omegagamma_plot_list = Float64[];

#         # Look at a single omega gamma value since each one spans all seeds
#         iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

#         for jdx in eachindex(iloop_omega_gamma_list)

#             # Get the idex for the current omega gamma value
#             matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
#             jloop_data_frame = iloop_data_frame[matching_jdx,:]

#             # get the mean over all seeds
#             jvalue_mean_wavespeed = mean(jloop_data_frame.wavespeed_x)

#             # Append values only if they are positive
#             if jvalue_mean_wavespeed >= 0
#                 push!(loop_mean_wavespeed_list, jvalue_mean_wavespeed)
#                 push!(loop_omegagamma_plot_list, iloop_omega_gamma_list[jdx])
#             end
#         end
#         scatterlines!(ax, loop_omegagamma_plot_list, loop_mean_wavespeed_list, color=marker_color, label=L"\hat{P} = %$(iloop_pressure_value), \hat{\gamma}=%$(gamma_value)")

#     end
#     axislegend(position = :lt)
#     fig
# end


function plot_ellipse_low_pressure(data_frame, gamma_values, pressure_value)

    # Define parameters to plot
    pressure_list = sort(unique(data_frame.input_pressure))
    closest_pressure_match_index = argmin(abs.(pressure_list .- pressure_value))
    plot_pressure = pressure_list[closest_pressure_match_index]
    gamma_list = sort(unique(data_frame.gamma))
    plot_gamma = [gamma_list[argmin(abs.(gamma_list .- element))] for element in gamma_values]


    # Filter the table to only those data
    matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
    matching_pressure_index = in.(data_frame.input_pressure, Ref(plot_pressure))
    combined_index = matching_gamma_index .& matching_pressure_index
    filtered_data_frame = data_frame[combined_index, :]

    # Limit range to data (deprecated)
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]

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
    set(ax_attenuation, 'XTickLabel', []);

    % Axes for Rotation Angle
    ax_rotation = nexttile;
    hold(ax_rotation, 'on');
    % xlabel(ax_rotation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel(ax_rotation, '\$ \\overline{\\left|\\theta\\right|} \$', "FontSize", 20, "Interpreter", "latex");
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
    normalized_variable = (log.(plot_gamma) .- minimum(log.(plot_gamma))) ./ (maximum(log.(plot_gamma)) .- minimum(log.(plot_gamma)))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_gamma)

        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # For idx, only show current gamma data
        iloop_gamma_value = plot_gamma[idx]
        iloop_gamma_index = in.(filtered_data_frame.gamma, Ref(iloop_gamma_value))
        iloop_combined_index = iloop_gamma_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_rotation_angles = Float64[];
        loop_mean_attenuation_list = Float64[];
        loop_mean_aspect_ratio_list = Float64[];

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
        valid_indices = iloop_omega_gamma_list .<= iloop_gamma_value.*2
        iloop_omega_gamma_list = iloop_omega_gamma_list[valid_indices]
        loop_mean_aspect_ratio_list = loop_mean_aspect_ratio_list[valid_indices]
        loop_mean_rotation_angles = loop_mean_rotation_angles[valid_indices]
        loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]
        
        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$ \\hat{\\gamma} = %.3f, \\hat{P} = %.3f\$", iloop_gamma_value, pressure_value)

        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(iloop_omega_gamma_list);
        mean_aspect_ratio = $(loop_mean_aspect_ratio_list);
        mean_rotation_angles = $(loop_mean_rotation_angles);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_gamma_value = $(iloop_gamma_value);
        plot_pressure = $(plot_pressure);
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

end