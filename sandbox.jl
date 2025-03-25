include("src/GranMA.jl")

using .GranMA
using DataFrames
using Debugger
using MATLAB
using Printf
using Statistics
using StatsBase
using Plots
using MAT
using Polynomials
using LinearAlgebra

simulation_data = load_data("out/processed/2d_bi_K100_W5.jld2")
simulation_data = load_data("out/processed/2d_K100_80kby40.jld2") # large packing with low freqs

# data_gaus = loadGausData("out/processed/gausGetAmps_noBacktrackV3.jld2")
data_gaus = loadGausData("out/processed/2d_gaus_withDamping.jld2")
# data_gaus = filter(x -> x.omega >= 0.03, data_gaus) # This is used for the presentation of the program review 
data_gaus = filter(x -> x.omega > .02, data_gaus) # This is for the wavespeed plot
data_gaus = filter(x -> x.pressure >= .001, data_gaus) 

function threeD()
    simulation_data = loadData3d("out/processed/3d_36900by7_V2.jld2") # 5by7by7
    simulation_data = loadData3d("out/processed/3d_80Kby15_V3.jld2") # 15by15 tiles out to 300 long channel

    plot_ωγ_attenuation_2d(simulation_data, .5, 1.2)

    # Mean Field Plots
    transverse_axis = "y"
    filtered_data = FilterData3d(simulation_data, .1, :pressure, .1, :omega, .5, :gamma, 1, :seed)
    getMeanField3d(filtered_data, transverse_axis)

    filtered_data = FilterData3d(simulation_data, .1, :pressure, .001, :omega, .5, :gamma, 1, :seed)
    
    # ----------------------- Final Equation Plots ------------------------
    plotAmp(filtered_data, transverse_axis="y")
    gamma_values = [ .001, .05, .1]
    plotStitchPhaseScatter3d(simulation_data, gamma_values, shear=true) 
    plotStitchAmpRatio3d(simulation_data, gamma_values)
end

function shearPaperPlots()
    # Shear data
    simulation_data = load_data("out/processed/2d_K100_shear_1000by5_all.jld2")
    simulation_data = load_data("out/processed/2d_shear_80kby40.jld2")

    # high omega
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .05, :gamma, 1, :seed)
    plotAmp(data, shear=true)
    plotPhase(data, shear=true) # phase plot for high pressure, low gamma
    data = FilterData(simulation_data, .001, :pressure, .02, :omega, .02, :gamma, 1, :seed)
    plotAmp(data, shear=true)
    plotPhase(data, shear=true) # phase plot for low pressure, low gamma

    # low omega
    data = FilterData(simulation_data, .1, :pressure, .001, :omega, .5, :gamma, 1, :seed)
    plotPhase(data, shear=true) # phase plot for high pressure, low gamma
    data = FilterData(simulation_data, .001, :pressure, .001, :omega, .5, :gamma, 1, :seed)
    plotPhase(data, shear=true) # phase plot for low pressure, low gamma
    
    # Mean field plots combined
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .1, :gamma, 1, :seed) # low pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data, shear = true) # save this as "f1.fig"
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .1, :gamma, 1, :seed) # high pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data, shear = true) # save this as "f1.fig"
 
    # mat"""addpath('src/matlab_functions'); combinePlotsTiled("f1.fig", "f2.fig", [0,200], [1E-2, 1])""" # Need to normalize the mean field
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("f1.fig", "f2.fig", [0,50], [1E-2, 1])""" # Need to normalize the mean field
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/f1_shear40wide.fig", "figures/f2_shear40wide.fig", [0,50], [1E-2, 1])""" # Need to normalize the mean field
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/f3_shear40wide.fig", "figures/f4_shear40wide.fig", [0,300], [1E-2, 1])""" # Need to normalize the mean fieldh
    # Single simulation
    mat"""addpath('src/matlab_functions'); sim2dShear(100, 1, 1, 1, 5000, 0.001, 5, 1, "in/2d_5wide_1000long/", "out/junk_yard")"""

    # Theory Plots
    gamma_values = [ .05, .1, .5, 1]
    gamma_values = [ .001, .05, .1]
    plotStitchPhaseScatter(simulation_data, gamma_values, shear=true) 
    plotStitchAmpRatio(simulation_data, gamma_values, shear=true)
    # plotStitchAmpPhase(simulation_data, gamma_values)
end

function paperPlots()
    sim2d(100, 1, 5, 1, 5000, .1, 5, 1)  # Single particle oscilaltion

    simulation_data = load_data("out/processed/2d_bi_K100_W5.jld2")
    # high omega
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .5, :gamma, 1, :seed)
    plotAmp(data) # amplitude plot for high pressure, low gamma
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .5, :gamma, 1, :seed)
    plotAmp(data) # amplitude plot for low pressure, low gamma
    
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .5, :gamma, 1, :seed)
    plotPhase(data) # phase plot for high pressure, low gamma
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .5, :gamma, 1, :seed)
    plotPhase(data) # phase plot for low pressure, low gamma

    # low omega
    data = FilterData(simulation_data, .1, :pressure, .01, :omega, .5, :gamma, 1, :seed)
    plotAmp(data) # amplitude plot for high pressure, low gamma
    data = FilterData(simulation_data, .001, :pressure, .01, :omega, .5, :gamma, 1, :seed)
    plotAmp(data) # amplitude plot for low pressure, low gamma

    data = FilterData(simulation_data, .1, :pressure, .01, :omega, .5, :gamma, 1, :seed)
    plotPhase(data) # phase plot for high pressure, low gamma
    data = FilterData(simulation_data, .001, :pressure, .01, :omega, .5, :gamma, 1, :seed)
    plotPhase(data) # phase plot for high pressure, low gamma

    # Tiled Plots
    # high omega
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .5, :gamma, 1, :seed) # High Pressure
    plotAmpTiled(data)
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .5, :gamma, 1, :seed) # low pressure
    plotAmpTiled(data)
    #low omega
    data = FilterData(simulation_data, .1, :pressure, .01, :omega, .5, :gamma, 1, :seed) # high pressure
    plotAmpTiled(data)
    data = FilterData(simulation_data, .001, :pressure, .01, :omega, .5, :gamma, 1, :seed) # low pressuer
    plotAmpTiled(data)
 
    # Amp and Phase Ratio Plots
    plotAmpRatio(simulation_data, .5) # save as fig1.fig
    plotPhaseRatio(simulation_data, .5) # save as fig2.fig
    mat"addpath('src/matlab_functions'); combinePlotsTiled('fig1_0.01.fig', 'fig2_0.01.fig')" # better to use the matlab version...
    mat"addpath('src/matlab_functions'); combinePlotsTiled('fig1_0.05.fig', 'fig2_0.05.fig')" # better to use the matlab version...
    mat"addpath('src/matlab_functions'); combinePlotsTiled('fig1_0.1.fig', 'fig2_0.1.fig')" # better to use the matlab version...
    mat"addpath('src/matlab_functions'); combinePlotsTiled('fig1_0.5.fig', 'fig2_0.5.fig')" # better to use the matlab version...
    mat"addpath('src/matlab_functions'); combinePlotsTiled('fig1.fig', 'fig2.fig')" # better to use the matlab version...

    plot_ωγ_attenuation_2d(simulation_data, .5, 1.2)
    data_gaus = loadGausData("out/processed/gausGetAmps_noBacktrackV3.jld2")
    plotGausAttenuation2d(data_gaus)

    # ellipse Plots
    sim2d(100, 1, 5, 1, 5000, .1 ,5,  1) # high pressure, 
    sim2d(100, 1, 5, 1, 5000, .001 ,5,  1)
    plotEllipseAttenuation2d(simulation_data, .5)
    
    # Derek's Mean field theory for single sim ---------------------------------------
    # These two below are have phase aligned and are used int the paper. Don't mess with them.
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .5, :gamma, 1, :seed) # low pressure,
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data)

    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .5, :gamma, 1, :seed) # high pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data)

    # Mean field plots combined
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .1, :gamma, 1, :seed) # low pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data) # save this as "f1.fig"
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .1, :gamma, 1, :seed) # high pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data) # save fiths as "f2.fig" 
 
    # mat"""addpath('src/matlab_functions'); combinePlotsTiled("f1.fig", "f2.fig", [0,200], [1E-2, 1])""" # Need to normalize the mean field
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/f1_compression_40wide_lowP.fig", "figures/f2_compression_40wide_lowP.fig", [0,200], [1E-2, 1])""" # Need to normalize the mean field
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/f1_compression_40wide_highP.fig", "figures/f2_compression_40wide_highP.fig", [0,200], [1E-2, 1])""" # Need to normalize the mean fieldh
    # 5 wide mean field
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/f1_compression_5wide_highP.fig", "figures/f2_compression_5wide_highP.fig", [0,200], [1E-2, 1])""" # Need to normalize the mean field
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/f1_compression_5wide_lowP.fig", "figures/f2_compression_5wide_lowP.fig", [0,200], [1E-2, 1])""" # Need to normalize the mean fieldh

    # Mean field but for shear  simulation_data = load_data("out/processed/2d_K100_shear_1000by5.jld2")
    data = FilterData(simulation_data, .1, :pressure, .001, :omega, .5, :gamma, 1, :seed) # high pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data, shear = true)

    plotAmp(data) # amplitude plot for low pressure, low gamma

    # mean field over all simulations
    plotAmpRatioMeanField(simulation_data, .1)  # save the figure as fig3.fig
    plotPhaseRatioMeanField(simulation_data, .1) # save this one too as fig4.fig
    mat"addpath('src/matlab_functions'); combinePlots('fig1.fig', 'fig3.fig')"
    mat"addpath('src/matlab_functions'); combinePlots('fig2.fig', 'fig4.fig')"
    # ----------------------------------------------------------------------------------

    # Stiched Equation Plots (testing the theory at begining of paper)
    gamma_values = [ .05, .1, .5, 1]
    plotStitchPhaseScatter(simulation_data, gamma_values) 
    plotStitchAmpRatio(simulation_data, gamma_values)
    plotStitchAmpPhase(simulation_data, gamma_values)

    gamma_values = [.08, .1,.3, .5, .6, .8, 1]
    plotStitchAttenuation(simulation_data, gamma_values, 1.2)
    plotStitchAttenuation(simulation_data, gamma_values, 1.2)

end

function width_effect_tile_test()
    # Smaller tiles (20 by 20)
    gamma_value = .6
    simulation_data = load_data("out/processed/2d_K100_tileTest_20by20.jld2")
    plot_ωγ_attenuation_2d(simulation_data, gamma_value, 1.2)

    # Medium tiles (40 by 20)
    simulation_data = load_data("out/processed/2d_K100_tileTest_40by20.jld2")
    plot_ωγ_attenuation_2d(simulation_data, gamma_value, 1.2)

    # Larger tiles (100 by 20)
    simulation_data = load_data("out/processed/2d_K100_tileTest_100by20.jld2")
    plot_ωγ_attenuation_2d(simulation_data, gamma_value, 1.2)

    # For the tile-test plots in the paper, save the 20by20 and 100by20 plots and use the following command
    # make sure to go into plotting.jl and comment out the loglog line in plot_ωγ_attenuation_2d to get rid of 1-D theory line
    mat"""addpath('src/matlab_functions'); combinePlots("figures/tileTest_20by20.fig", "figures/tileTest_100by20.fig")""" 
end



# Ellipse

function plot_ellipse_width_effect(ω_value, γ_value, pressure_value)
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
        filtered_data = FilterData(simulation_data, pressure_value, :pressure)
        probabilities_asp , probabilities_rot, plot_bins_asp, plot_bins_rot = plot_ellipse_pdf(ω_value, γ_value, plot=false, simulation_data = filtered_data) 

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

function plotEllipseFlattness(γ_value; simulation_data = simulation_data)
    filtered_data = FilterData(simulation_data, γ_value, :gamma)
    pressure_list = sort(unique([entry.pressure for entry in filtered_data])) # goes through each entry of simulation_data and get the P value at that entry
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))
    mat"""
    figure_aspect_ratio = figure;
    xlabel('\$ \\hat{\\omega} \\hat{\\gamma} \$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\textrm{Flattness} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    hold on
    """
    for pressure_value in pressure_list
        # Only look at data for current pressure value
        pressureLoop_filtered_data = FilterData(filtered_data, pressure_value, :pressure)
    
        # Assign a color
        idx = findfirst(idx -> idx == pressure_value, pressure_list)
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]
    
        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in pressureLoop_filtered_data]))
        slope_list = Float64[]
        
        for omega_gamma_value in matching_omega_gamma_list
    
            OmegaGamma_filtered_data = FilterData(pressureLoop_filtered_data, omega_gamma_value, :omega_gamma)
            ω_value = OmegaGamma_filtered_data[1].omega
    
            # Look at a single omega gamma value since each one spans all seeds
            seed_list = sort(unique([entry.seed for entry in OmegaGamma_filtered_data]))
    
            slope_asp = Float64[]
            slope_rot = Float64[]
    
            for seed_value in seed_list
                seedLoop_filtered_data = FilterData(OmegaGamma_filtered_data, seed_value, :seed)
                probabilities_asp , probabilities_rot, plot_bins_asp, plot_bins_rot = plot_ellipse_pdf(ω_value, γ_value; plot=false, simulation_data=seedLoop_filtered_data)
    
            # Flatten the vectors
            probabilities_asp = vec(probabilities_asp)
            probabilities_rot = vec(probabilities_rot)
            plot_bins_asp = vec(plot_bins_asp)
            plot_bins_rot = vec(plot_bins_rot)

            # Align lengths by trimming to the minimum length
            min_length_asp = min(length(probabilities_asp), length(plot_bins_asp))
            min_length_rot = min(length(probabilities_rot), length(plot_bins_rot))

            probabilities_asp = probabilities_asp[1:min_length_asp]
            plot_bins_asp = plot_bins_asp[1:min_length_asp]

            probabilities_rot = probabilities_rot[1:min_length_rot]
            plot_bins_rot = plot_bins_rot[1:min_length_rot]



                slope_asp_value = calculate_slope(plot_bins_asp, probabilities_asp)
                slope_asp_value = isnan(slope_asp_value) ? NaN : (isinf(slope_asp_value) ? NaN : slope_asp_value)
    
                # Calculate slope for rotational probability (probabilities_rot vs plot_bins_rot)
                slope_rot_value = calculate_slope(plot_bins_rot, probabilities_rot)
                slope_rot_value = isnan(slope_rot_value) ? NaN : (isinf(slope_rot_value) ? NaN : slope_rot_value)
    
                # Store the slopes for both asp and rot
                slope_asp = push!(slope_asp, slope_asp_value)
                slope_rot = push!(slope_rot, slope_rot_value)
                println("$slope_asp")
            end
            
            # Average slope across seeds and push to slope_list
            slope_list = push!(slope_list,  ( mean(slope_asp) * mean(slope_rot)))
        end
        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$ \\hat{P} = %.3f, \\hat{\\gamma} = %.3f  \$", pressure_value, γ_value)

        mat"""
        y = $(slope_list);
        x = $(matching_omega_gamma_list);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);
        
        figure(figure_aspect_ratio);
        plot(x, y, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);

        """

    end
    mat"""
    legend('show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    """
end

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
    pressure_list = sort(unique([entry.pressure for entry in matching_ωγ_data])) # goes through each entry of simulation_data and get the P value at that entry
    
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
        xlabel('\$ \\overline{ \\left| \\theta \\right| } \$', "FontSize", 20, "Interpreter", "latex");
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


# Attenuation Plots

function plotOmegaAttenuation2d(simulation_data; plot=true)  
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
            loop_mean_alphaoveromega = mean(filter(x -> x > 0, [entry.attenuation_x for entry in matching_omega_data])) ./ omega_value
            loop_std_alphaoveromega =  std(entry.attenuation_x for entry in matching_omega_data) ./ omega_value
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

function plot_γ_attenuation_2d(ω_value; plot=true, simulation_data=simulation_data)  # Need to fix the legend
    
    filtered_data = FilterData(simulation_data, ω_value, :omega)


    pressure_list = sort(unique([entry.pressure for entry in filtered_data])) # goes through each entry of simulation_data and get the P value at that entry

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);

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
        pressureLoop_filtered_data = FilterData(filtered_data, pressure_value, :pressure )
        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in pressureLoop_filtered_data]))

        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current pressure value
            omegaGammaLoop_filtered_data = FilterData(pressureLoop_filtered_data, omega_gamma_value, :omega_gamma)

            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in omegaGammaLoop_filtered_data)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        # valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        # matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        # loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]

        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{\\omega} = %.3f, \\hat{P} = %.3f \$", ω_value, pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(matching_omega_gamma_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            plot_gamma = $(ω_value);
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
end

# migrated to plotting.jl
# function plot_ωγ_attenuation_2d(gamma_value, mean_diameter; plot=true, simulation_data=simulation_data)  # Need to fix the legend
#     # Initialize outputs
#     matching_omega_gamma_list = []
#     loop_mean_attenuation_list = []

#     # filter the data based on those that are close to gamma_value
#     closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
#     closest_gamma_value = simulation_data[closest_gamma_index].gamma
#     matching_gamma_data = filter(entry -> entry.gamma == closest_gamma_value, simulation_data)
#     plot_gamma = gamma_value

#     # Get a list of unique input pressures
#     pressure_list = sort(unique([entry.pressure for entry in matching_gamma_data])) # goes through each entry of simulation_data and get the P value at that entry
#     plot_pressure = pressure_list

#     # Define the plot limits to match the 1D theory plot curves
#     theory_x = collect(3E-4:1E-5:3)
#     theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
#     # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
#     # upper_limit_line_y = [1E-5; 1]
#     # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
#     # lower_limit_line_y = [1E-5; 1]
#     # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
#     # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

#     if plot
#         # Intialized the plots to loop over
#         mat"""
#         figure_attenuation = figure;
#         loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
#         hold on;
#         xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
#         ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
#         set(gca, 'XScale', 'log');
#         set(get(gca, 'ylabel'), 'rotation', 0);
#         grid on;
#         box on;
#         """
#     end

#     # Normalize the gamma values
#     normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

#     # Create a line for each gamma value across all pressure_list
#     for pressure_value in pressure_list

#         # Assign a color
#         idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
#         marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

#         # Only look at data for current pressure value
#         matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

#         # Initizalized vectors for just this pressure
#         loop_mean_attenuation_list = Float64[];

#         # Look at a single omega gamma value since each one spans all seeds
#         matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))

#         for omega_gamma_value in matching_omega_gamma_list

#             # Only look at data for current pressure value
#             matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

#             # Get the mean over all seeds
#             loop_mean_alphaoveromega = mean_diameter .* mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)

#             # Append values
#             push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
#         end

#         # Filter data to include only points where omega_gamma <= gamma_value
#         valid_indices = matching_omega_gamma_list .<= gamma_value.*2
#         matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
#         loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]

#         # This is needed because MATLAB.jl has a hard time escaping \'s
#         legend_label = @sprintf("\$ \\hat{P} = %.3f, \\hat{\\gamma} = %.3f\$", pressure_value, gamma_value)

#         if plot
#             # Transfer data to MATLAB
#             mat"""
#             omega_gamma = $(matching_omega_gamma_list);
#             mean_attenuation_x = $(loop_mean_attenuation_list);
#             iloop_pressure_value = $(pressure_value);
#             plot_gamma = $(plot_gamma);
#             marker_color= $(marker_color);
#             legend_label = $(legend_label);
#             figure(figure_attenuation);
#             set(gca, 'Yscale', 'log');
#             plot(omega_gamma, mean_attenuation_x, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
#             """
#         end
#     end
#     if plot
#         # Add legends to the plots
#         mat"""
#         % Add legends to the MATLAB plots
#         figure(figure_attenuation);
#         legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
#         """
#     end
#     return matching_omega_gamma_list, loop_mean_attenuation_list
# end

function plot_γ_attenuation_P_2d(pressure_value, ω_values::Vector{Float64}; plot=true, simulation_data=simulation_data)

    if plot

        theory_x = collect(3E-4:1E-5:3)
        theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
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
    
    closest_ω_indices = [argmin(abs.([idx.omega for idx in simulation_data] .- ω)) for ω in ω_values]
    closest_ω_list = [simulation_data[idx].omega for idx in closest_ω_indices]
    normalized_variable = (log.(closest_ω_list) .- minimum(log.(closest_ω_list))) ./ (maximum(log.(closest_ω_list)) .- minimum(log.(closest_ω_list)))


    for ω_value in closest_ω_list
        idx = findfirst(idx -> idx ==ω_value, closest_ω_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]
        filtered_data = FilterData(simulation_data, ω_value, :omega, pressure_value, :pressure)
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in filtered_data]))
        loop_mean_attenuation_list = Float64[];

        for omega_gamma_value in matching_omega_gamma_list
            # Only look at data for current omega_gamma_value 
            omegaGammaLoop_filtered_data = FilterData(filtered_data, omega_gamma_value, :omega_gamma)

            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in omegaGammaLoop_filtered_data)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
        end
        legend_label = @sprintf("\$ \\hat{\\omega} = %.3f, \\hat{P} = %.3f \$", ω_value, pressure_value)

        if plot
            mat"""
            omega_gamma = $(matching_omega_gamma_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            plot_gamma = $(ω_value);
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
end

# migrated to plotting.jl
# function plot_ωγ_wavespeed_2d(gamma_value) # Need to fix lgend

#     # filter the data based on those that are close to gamma_value
#     closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
#     closest_gamma_value = simulation_data[closest_gamma_index].gamma
#     matching_gamma_data = filter(entry -> entry.gamma == closest_gamma_value, simulation_data)
#     plot_gamma = gamma_value

#     # Get a list of unique input pressures
#     pressure_list = unique([entry.pressure for entry in matching_gamma_data]) # goes through each entry of simulation_data and get the P value at that entry
#     plot_pressure = pressure_list

#     # Define the plot limits to match the 1D theory plot curves
#     theory_x = collect(3E-4:1E-5:3)
#     theory_y = 1 ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(0.5) ./ (1 .+ theory_x.^2);
#     # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
#     # upper_limit_line_y = [1E-5; 1]
#     # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
#     # lower_limit_line_y = [1E-5; 1]
#     # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
#     # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

#     # Intialized the plots to loop over
#     mat"""
#     figure_wavespeed = figure;
#     loglog($(theory_x), 1./$(theory_y), 'k', 'DisplayName', '1-D Theory'), hold on
#     xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
#     ylabel('\$\\hat{c} \$', "FontSize", 20, "Interpreter", "latex");
#     set(gca, 'XScale', 'log');
#     set(get(gca, 'ylabel'), 'rotation', 0);
#     grid on;
#     box on;
#     """

#     # Normalize the gamma values
#     normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

#     # Create a line for each gamma value across all pressure_list
#     for pressure_value in pressure_list

#         # Assign a color
#         idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
#         marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

#         # Only look at data for current pressure value
#         matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

#         # Initizalized vectors for just this pressure
#         loop_mean_wavespeed_list = Float64[];

#         # Look at a single omega gamma value since each one spans all seeds
#         matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))

#         for omega_gamma_value in matching_omega_gamma_list

#             # Only look at data for current pressure value
#             matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
#             @bp
#             # Get the mean over all seeds
#             loop_mean_wavespeed = mean(entry.:wavespeed_x for entry in matching_omega_gamma_data)
#             @bp
#             # Append values
#             push!(loop_mean_wavespeed_list, loop_mean_wavespeed)
#         end

#         # Filter data to include only points where omega_gamma <= gamma_value
#         valid_indices = matching_omega_gamma_list .<= gamma_value.*2
#         matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
#         loop_mean_wavespeed_list = loop_mean_wavespeed_list[valid_indices]
        
#         mat"""
#         omega_gamma = $(matching_omega_gamma_list);
#         mean_wavespeed_x = $(loop_mean_wavespeed_list);
#         iloop_pressure_value = $(pressure_value);
#         plot_gamma = $(plot_gamma);
#         marker_color = $(marker_color);
#         pressure_label = sprintf('Wavespeed X = %.2f, Gamma = %.2f (Aspect Ratio)', $(pressure_value), $(plot_gamma));
#         % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
#         pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Wavespeed X)', $(pressure_value), $(plot_gamma));
        
#         figure(figure_wavespeed);
#         set(gca, 'Yscale', 'log');
#         plot(omega_gamma, mean_wavespeed_x, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
#         """
#     end

#     # Add legends to the plots
#     mat"""
#     % Add legends to the MATLAB plots
#     figure(figure_wavespeed);
#     legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
#     """
# end


function plot_attenuation_width_effect(γ_value, pressure_value)
    width_list = [10, 20 , 50]

    normalized_variable = (width_list .- minimum(width_list)) ./ (maximum(width_list) .- minimum(width_list))

    # Intialized the plots to loop over
    mat"""
    figure_attenuation = figure;
    xlabel('\$ \\hat{\\omega}\\hat{\\gamma} \$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log')
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    hold on

    """
    
    for width in width_list
        @bp
        # Assign a color
        idx = findfirst(idx -> idx ==width, width_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]
        
        simulation_data = load_data("out/processed/2d_bi_K100_W$(width).jld2")
        filtered_data = FilterData(simulation_data, pressure_value, :pressure)
        matching_omega_gamma_list, loop_mean_attenuation_list = plot_ωγ_attenuation_2d(γ_value, plot=false, simulation_data=filtered_data)
        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.3f, \\textrm{Width} = %.3f, \\hat{\\gamma} = %.3f \$", pressure_value,width, γ_value)

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

function bin_plot_energy(pressure_value, γ_value, ω_value, seed_value; plot=true, simulation_data=simulation_data)

        
    filtered_data = FilterData(simulation_data, pressure_value, :pressure, γ_value, :gamma, ω_value, :omega, seed_value, :seed)

    omega = ω_value
    gamma = γ_value
    
    # Get input for this simulation
    amp_vector_y = filtered_data[1].amplitude_vector_y
    amp_vector_x = filtered_data[1].amplitude_vector_x
    phase_vector_y = filtered_data[1].unwrapped_phase_vector_y 
    phase_vector_x = filtered_data[1].unwrapped_phase_vector_x 
    distance_vector_y = filtered_data[1].initial_distance_from_oscillation_output_y_fft 
    distance_vector_x = filtered_data[1].initial_distance_from_oscillation_output_x_fft 

    # Check if any of the vectors are empty
    if isempty(amp_vector_y) || isempty(phase_vector_y) || isempty(distance_vector_y)
        println("Warning: One or more input vectors are empty. Returning NaN for Q_ratio.")
        return NaN
    end

    # Using fit(Histogram) to divide distance_vector_y from 1 to max distance away from wall
    bins_y = StatsBase.fit(Histogram, distance_vector_y, 1:maximum(distance_vector_y)+1)
    bins_x = StatsBase.fit(Histogram, distance_vector_x, 1:maximum(distance_vector_x)+1)

    # Initialize vectors to store bin centers and energy losses for plotting
    bin_centers_y = Float64[]
    energy_losses_y = Float64[]
    bin_centers_x = Float64[]
    energy_losses_x = Float64[]

    # Assign a total energy loss for each bin
    for i in 1:length(bins_y.weights)
        indices = findall(bins_y.edges[1][i] .<= distance_vector_y .< bins_y.edges[1][i+1])
        bin_center = (bins_y.edges[1][i] + bins_y.edges[1][i+1]) / 2
        if length(indices) > 1
            amp_bin = amp_vector_y[indices]
            phase_bin = phase_vector_y[indices]
            energy_loss = total_energy_loss_in_bin(amp_bin, phase_bin, gamma, omega)
            if energy_loss > 0  # Only include positive energy losses
                push!(bin_centers_y, bin_center)
                push!(energy_losses_y, energy_loss)
            end
        end
    end

    for i in 1:length(bins_x.weights)
        indices = findall(bins_x.edges[1][i] .<= distance_vector_x .< bins_x.edges[1][i+1])
        bin_center = (bins_x.edges[1][i] + bins_x.edges[1][i+1]) / 2
        if length(indices) > 1
            amp_bin = amp_vector_x[indices]
            phase_bin = phase_vector_x[indices]
            energy_loss = total_energy_loss_in_bin(amp_bin, phase_bin, gamma, omega)
            if energy_loss > 0  # Only include positive energy losses
                push!(bin_centers_x, bin_center)
                push!(energy_losses_x, energy_loss)
            end
        end
    end

    # Perform a linear fit on the log-transformed energy losses
    if !isempty(bin_centers_y) && !isempty(energy_losses_y)
        log_energy_losses_y = log.(energy_losses_y)
        p_y = Polynomials.fit(bin_centers_y, log_energy_losses_y, 1)
        slope_y = coeffs(p_y)[2]  # Get the slope from the coefficients
        intercept_y = coeffs(p_y)[1]
        fitted_y = exp.(intercept_y .+ slope_y .* bin_centers_y)  # Convert back to linear space
    else
        intercept_y, slope_y, fitted_y = NaN, NaN, []
    end

    if !isempty(bin_centers_x) && !isempty(energy_losses_x)
        log_energy_losses_x = log.(energy_losses_x)
        p_x = Polynomials.fit(bin_centers_x, log_energy_losses_x, 1)
        slope_x = coeffs(p_x)[2]  # Get the slope from the coefficients
        intercept_x = coeffs(p_x)[1]
        fitted_x = exp.(intercept_x .+ slope_x .* bin_centers_x)  # Convert back to linear space
    else
        intercept_x, slope_x, fitted_x = NaN, NaN, []
    end

    if plot
        # Plotting the results using MATLAB

        title_label = @sprintf("\$\\hat{P} = %.4f, \\hat{\\gamma} = %.4f, \\hat{\\omega} = %.4f \$", pressure_value, γ_value, ω_value)
        mat"""
        title_label = $(title_label)
        set(groot, 'defaultTextInterpreter', 'latex');  % Set for text objects
        set(groot, 'defaultLegendInterpreter', 'latex');  % Set for legends
        set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  % Set for axes tick labels
        set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');  % Set for colorbar tick labels
        set(groot, 'defaultTextarrowshapeInterpreter', 'latex');  % Set for text arrows in annotations
        
        % Plot for y-direction
        scatter($(bin_centers_y), $(energy_losses_y), 'DisplayName', 'y')
        hold on
        if length($(fitted_y)) > 0
            plot($(bin_centers_y), $(fitted_y), 'DisplayName', sprintf('Fit y (slope = %.3f)', $(slope_y)))
        end

        % Plot for x-direction
        scatter($(bin_centers_x), $(energy_losses_x), 'DisplayName', 'x')
        if length($(fitted_x)) > 0
            plot($(bin_centers_x), $(fitted_x), 'DisplayName', sprintf('Fit x (slope = %.3f)', $(slope_x)))
        end

        % Set logarithmic scale, labels, and other plot settings
        set(gca, 'yscale', 'log')
        xlabel(' \$ x- \$ Distance from Oscillation', 'FontSize', 15)
        ylabel(' \$ \\left< E \\right> \$', 'Rotation', 0, 'FontSize', 15)
        grid on
        box on
        legend()
        title($title_label)
        """
    end
    # Q_ratio = exp(intercept_y) #/ ω_value
    Q_ratio =  - slope_y 
    return Q_ratio 
end

function total_energy_loss_in_bin(amps, phases, gamma, omega)
    N = length(amps)
    total_loss = 0.0
    for j in 1:(N-1)
        total_loss += energy_loss_time_avg(amps[j], phases[j], amps[j+1], phases[j+1], gamma, omega)
    end
    return total_loss
end

# This function assumes a spring constant K=100
function energy_loss_time_avg(amp_j, phi_j, amp_k, phi_k, gamma, omega)
    K = 100
    m = 1
    # K*omega/2 * amp_j * amp_k * sin(phi_j - phi_k) -  (gamma * omega^2 )*((amp_j^2 + amp_k^2)/2 - amp_j * amp_k * cos(phi_j - phi_k))
    return abs(omega * sqrt(K^3 / m) * ( .5 * amp_j * amp_k * sin(phi_j - phi_k) - (gamma * omega ) * ((amp_j^2 + amp_k^2)/2 - amp_j * amp_k * cos(phi_j - phi_k))))
end

function plot_energy(γ_value) 

    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
    plot_gamma = γ_value
    gamma_value = γ_value

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry

    # Limit range to data
    upper_limit_line_x = [1*γ_value; 1*γ_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*γ_value; .1*γ_value]
    lower_limit_line_y = [1E-5; 1]

    # Start MATLAB session
    mat"""
    figure_main = figure;
    tiled_main = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'none'); % 2 rows, 1 column

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

    % Axes for Energy
    ax_energy = nexttile;
    hold(ax_energy, 'on');
    % xlabel(ax_energy, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel(ax_energy, '\$ \\frac{ \\overline{E} }{\\hat{\\omega}}  \$', "FontSize", 20, "Interpreter", "latex");
    set(ax_energy, 'XScale', 'log');
    set(ax_energy, 'YScale', 'log')
    set(get(ax_energy, 'ylabel'), 'rotation', 0);
    grid(ax_energy, 'on');
    box(ax_energy, 'on');
    set(ax_energy, 'XTickLabel', []);
    """

    # get a range for plotting color from 0 to 1
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_E_list = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
        
        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current omega_gamma value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds       
            jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
            E_ratio_list = Float64[]
            seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))
            for k_seed in seed_list
                k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                k_seed_omega = k_seed_data[1].omega
                k_E_ratio = bin_plot_energy(pressure_value, γ_value, k_seed_omega, k_seed; plot=false, simulation_data=matching_omega_gamma_data)
                k_E_ratio = k_E_ratio / k_seed_omega
                push!(E_ratio_list, k_E_ratio)
            end
            j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
            push!(loop_mean_E_list, j_E_ratio)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        # valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        # matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        # loop_mean_E_list = loop_mean_E_list[valid_indices]
        # loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]
        @bp
        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$\\hat{P} = %.4f, \\hat{\\gamma} = %.2f\$", pressure_value, plot_gamma)

        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        loop_mean_E_list = $(loop_mean_E_list);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);

        % Plot Attenuation
        loglog(ax_attenuation, omega_gamma, mean_attenuation_x, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
        % Plot Aspect Ratio
        plot(ax_energy, omega_gamma, loop_mean_E_list, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        """
    end

    # Add legends to the plots
    mat"""
    legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    %legend(ax_energy, 'show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end




# Helper Functions

# function wrapped_distance(y1, y2)
#     direct_dist = abs(y1 - y2)
#     wrapped_dist = 2π - direct_dist
#     return min(direct_dist, wrapped_dist)
# end

# function mean_nearest_neighbor_distance(x_values, y_values)
#     points = hcat(x_values, y_values)
#     n = size(points, 1)
#     distances = zeros(n)
    
#     for i in 1:n
#         dist_to_others = zeros(n)
        
#         for j in 1:n
#             if i != j
#                 # Calculate Euclidean distance for x and wrapped distance for y
#                 dist_x = abs(x_values[i] - x_values[j])
#                 dist_y = wrapped_distance(y_values[i], y_values[j])
#                 dist_to_others[j] = sqrt(dist_x^2 + dist_y^2)
#             else
#                 dist_to_others[j] = Inf  # Exclude the point itself
#             end
#         end
        
#         # Find the minimum distance to the nearest neighbor
#         distances[i] = minimum(dist_to_others)
#     end
    
#     # Return the mean nearest neighbor distance
#     return mean(distances)
# end

function calculate_slope(x, y)
    p = Polynomials.fit(x, y, 1)  # Fit a 1st-degree polynomial (linear regression)
    return coeffs(p)[2]  # The slope is the coefficient of x
end