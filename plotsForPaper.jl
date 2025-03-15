include("src/GranMA.jl")

using .GranMA
using MATLAB

# 2D Compression Plots

function compressionMeanField2D()

    # 5 wide --------------------------------------------
    simulation_data = load_data("out/processed/2d_bi_K100_W5.jld2")

    # low pressure
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .1, :gamma, 1, :seed) # low pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data) # save this as "f1.fig"
    mat"saveas(figure(1), 'figures/meanField_compression_2d_5wide_lowP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_compression_2d_5wide_lowP_phase.fig')"
    mat"close all"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_compression_2d_5wide_lowP_phase.fig", "figures/meanField_compression_2d_5wide_lowP_amp.fig", [0,200], [1E-2, 1])""" 
    mat"saveas(figure(1), 'figures/meanField_compression_2d_5wide_lowP.eps')"
    mat"close all"

    # high pressure
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .1, :gamma, 1, :seed)
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data)
    mat"saveas(figure(1), 'figures/meanField_compression_2d_5wide_highP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_compression_2d_5wide_highP_phase.fig')"
    mat"close all"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_compression_2d_5wide_highP_phase.fig", "figures/meanField_compression_2d_5wide_highP_amp.fig", [0,200], [1E-2, 1])""" 
    mat"saveas(figure(1), 'figures/meanField_compression_2d_5wide.eps')"
    mat"close all"

    # 40 wide --------------------------------------------
    simulation_data = load_data("out/processed/2d_K100_80kby40.jld2") # large packing with low freqs

    # low pressure
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .1, :gamma, 1, :seed) # low pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data) # save this as "f1.fig"
    mat"saveas(figure(1), 'figures/meanField_compression_2d_40wide_lowP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_compression_2d_40wide_lowP_phase.fig')"
    mat"close all"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_compression_2d_40wide_lowP_phase.fig", "figures/meanField_compression_2d_40wide_lowP_amp.fig", [0,200], [1E-2, 1])""" 
    mat"saveas(figure(1), 'figures/meanField_compression_2d_40wide_lowP.eps')"
    mat"close all"

    # high pressure
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .1, :gamma, 1, :seed)
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data)
    mat"saveas(figure(1), 'figures/meanField_compression_2d_40wide_highP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_compression_2d_40wide_highP_phase.fig')"
    mat"close all"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_compression_2d_40wide_highP_phase.fig", "figures/meanField_compression_2d_40wide_highP_amp.fig", [0,200], [1E-2, 1])""" 
    mat"saveas(figure(1), 'figures/meanField_compression_2d_40wide.eps')"
    mat"close all"
end


# 2d Shear Plots

function shearMeanField2D()

    simulation_data = load_data("out/processed/2d_shear_80kby40.jld2")

    # low pressure
    data = FilterData(simulation_data, .001, :pressure, .1, :omega, .1, :gamma, 1, :seed) # low pressure
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data, shear = true)
    mat"saveas(figure(1), 'figures/meanField_shear_2d_40wide_lowP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_shear_2d_40wide_lowP_phase.fig')"
    mat"close all"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_shear_2d_40wide_lowP_phase.fig", "figures/meanField_shear_2d_40wide_lowP_amp.fig", [0,50], [1E-2, 1])""" 
    mat"saveas(figure(1), 'figures/meanField_shear_2d_40wide_lowP.eps')"
    mat"close all"

    # high pressure
    data = FilterData(simulation_data, .1, :pressure, .1, :omega, .1, :gamma, 1, :seed)
    mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(data, shear = true)
    mat"saveas(figure(1), 'figures/meanField_shear_2d_40wide_highP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_shear_2d_40wide_highP_phase.fig')"
    mat"close all"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_shear_2d_40wide_highP_phase.fig", "figures/meanField_shear_2d_40wide_highP_amp.fig", [0,300], [1E-2, 1])""" 
    mat"saveas(figure(1), 'figures/meanField_shear_2d_40wide.eps')"
    mat"close all"
end