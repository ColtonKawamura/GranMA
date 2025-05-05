include("src/GranMA.jl")

using .GranMA
using MATLAB
mat"set(0, 'DefaultFigureVisible', 'on')"
mat"set(groot, 'defaultFigureCreateFcn', @(fig,~)addlistener(fig, 'ObjectBeingDestroyed', @(~,~)disp(['Figure destroyed: ' num2str(fig.Number)])))"
mat"addpath('src/matlab_functions')"
mat"addpath('src/matlab_functions')"


# ----------------------------------------------------------------
# ------------------- 3D  ----------------------------------------
# ----------------------------------------------------------------
simulation_data = loadData3d("out/processed/3d_80Kby15_V4_1.jld2") # 15by15 tiles with yz intial positions
saveData3dToMat(simulation_data, "out/processed/3d_80Kby15_V4_1.mat") # save data to matlab format

#  3D Compression Attenuation  
    gamma_value = .2
    plot_ωγ_attenuation_2d(simulation_data, gamma_value, 1.2)    
    plot_ωγ_wavespeed_2d(simulation_data, gamma_value) 
    mat"slopeLine('loglog' ,1/4, [.05,.5], .3)" # for gamma = .5
    mat"slopeLine('loglog' ,1/6, [.02,.2], .35)" # for gamma = .2

#  3D Sitched Attenuation  
    gamma_values = [.1, .2, .5]
    plotStitchAttenuation(simulation_data, gamma_values, 1.2)
    mat"slopeLine('loglog' ,.75, [.005,.1], .05)"
    mat"slopeLine('loglog' ,.17, [.005,.1], .25)"
    

#  3D Compression MeanField 
    # Low pressure
    # z-axis
    filtered_data = FilterData3d(simulation_data, .001, :pressure, .1, :omega, .1, :gamma, 1, :seed, [5,10], "y", [0,2], "z")
    filtered_data = FilterData3d(simulation_data, .001, :pressure, .1, :omega, .1, :gamma, 1, :seed);
    transverse_axis = "z";
    getMeanField3d(filtered_data, transverse_axis);
    mat"saveas(figure(1), 'figures/meanField_compression_3d_z_15wide_lowP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_compression_3d_z_15wide_lowP_phase.fig')"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_compression_3d_z_15wide_lowP_phase.fig", "figures/meanField_compression_3d_z_15wide_lowP_amp.fig", [0,80], [1E-2, 1])""" 

    # y-ax
    transverse_axis = "y"
    getMeanField3d(filtered_data, transverse_axis)
    mat"saveas(figure(1), 'figures/meanField_compression_3d_y_15wide_lowP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_compression_3d_y_15wide_lowP_phase.fig')"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_compression_3d_y_15wide_lowP_phase.fig", "figures/meanField_compression_3d_y_15wide_lowP_amp.fig", [0,80], [1E-2, 1])""" 
    # mat"saveas(figure(1), 'figures/meanField_compression_3d_y_15wide_lowP.eps')"

    # High pressure
    # z-axis
    filtered_data = FilterData3d(simulation_data, .1, :pressure, .1, :omega, .1, :gamma, 1, :seed)
    transverse_axis = "z"
    getMeanField3d(filtered_data, transverse_axis)
    mat"saveas(figure(1), 'figures/meanField_compression_3d_z_15wide_highP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_compression_3d_z_15wide_highP_phase.fig')"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_compression_3d_z_15wide_highP_phase.fig", "figures/meanField_compression_3d_z_15wide_highP_amp.fig", [0,150], [1E-2, 1])""" 
    mat"saveas(figure(1), 'figures/meanField_compression_3d_z_15wide_highP.eps')"

    # y-ax
    transverse_axis = "y"
    getMeanField3d(filtered_data, transverse_axis);
    mat"saveas(figure(1), 'figures/meanField_compression_3d_y_15wide_highP_amp.fig')"
    mat"saveas(figure(2), 'figures/meanField_compression_3d_y_15wide_highP_phase.fig')"
    mat"""addpath('src/matlab_functions'); combinePlotsTiledTwo("figures/meanField_compression_3d_y_15wide_highP_phase.fig", "figures/meanField_compression_3d_y_15wide_highP_amp.fig", [0,150], [1E-2, 1])""" 
    mat"saveas(figure(1), 'figures/meanField_compression_3d_y_15wide_highP.eps')"


# 3D Energy per Cycle Theory
    filtered_data = FilterData3d(simulation_data, [5,10], "y", [0,2], "z")
    simulation_data = nothing
    GC.gc()
    gamma_values = [ .05, .1, .5, 1]
    plotStitchPhaseScatter3d(filtered_data, gamma_values);
    mat"slopeLine('loglog' ,.833, [.08,1], .1)" # (type, slope, [xlower, xupper], yMean)
    plotStitchPhaseScatter3d(simulation_data, gamma_values);
    mat"slopeLine('loglog' ,.833, [.08,1], .1)" # (type, slope, [xlower, xupper], yMean)
    
    # plotStitchAmpRatio(simulation_data, gamma_values); neeed 3d version
    # plotStitchAmpPhase(simulation_data, gamma_values);
# ----------------------------------------------------------------
# ------------------- 2D  ----------------------------------------
# ----------------------------------------------------------------

## 2D Compression
simulation_data = load_data("out/processed/2d_K100_80kby40.jld2") # large packing with low freqs
saveData2dToMat(simulation_data, "out/processed/2d_K100_80kby40.mat") # save data to matlab format
simulation_data = load_data("out/processed/2d_bi_K100_W5.jld2")
saveData2dToMat(simulation_data, "out/processed/2d_bi_K100_W5.mat") # save data to matlab format

# Shear Conversion
simulation_data = load_data("out/processed/2d_shear_80kby40.jld2")
saveData2dToMat(simulation_data, "out/processed/2d_shear_80kby40.mat") # save data to matlab format

# Compression Attenuation and Wavespeed
    gamma_value = .5
    plot_ωγ_attenuation_2d(simulation_data, gamma_value, 1.2)
    plot_ωγ_wavespeed_2d(simulation_data, gamma_value) 

# 2D Compression MeanField 

    simulation_data = load_data("out/processed/2d_bi_K100_W5.jld2") # 5 wide packing with low freqs

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

# 2D Compression Mode Density
    mat"""plotDampedModeDensityPDF("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N400_K100_M1.mat", [.2, .01, .001], [.1])"""

# 2D Stiched Attenution
    gamma_values = [.08, .1,.3, .5, .6, .8, 1]
    plotStitchAttenuation(simulation_data, gamma_values, 1.2)

# 2D Energy per Cycle Theory
    gamma_values = [ .05, .1, .5, 1]
    plotStitchPhaseScatter(simulation_data, gamma_values);
    plotStitchAmpRatio(simulation_data, gamma_values);
    plotStitchAmpPhase(simulation_data, gamma_values);


## 2D Shear
simulation_data = load_data("out/processed/2d_shear_80kby40.jld2")
# 2D Shear Mean Field
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


# ----------------------------------------------------------------
# ------------------- Packings -----------------------------------
# ----------------------------------------------------------------


# ------------------- 2d Poly Packings -------------------
# Pack Poly Tiles
    mat"packPoly2dRepXY(400,100, 1, 1, 1, .1, 20, 1, false, 1, 1, false, 'in/2d_poly_20by20/tiles/')"
# Repeat Poly Packings
    mat"pack2dRepeatTilePoly(400, 100, .1, 20, 1, 100, 2, false, 'in/2d_poly_20by20/tiles/', 'in/2d_poly_20by20/')"


# ------------------- Packing Visulaization -------------------
# Visulize Tile (not good for big tiles)
    mat"visualizeTile('in/2d_poly_20by20/tiles/2D_poly_N400_P0.001_Width20_Seed1.mat')"
# Visualize Full Packing
    mat"visualizeFullPacking('in/2d_poly_20by20/2D_N80000_P0.001_Width40_Seed1.mat', 20, 20)"