include("src/GranMA.jl")

using .GranMA
using MATLAB


# ----------------------------------------------------------------
# ------------------- 3D  -----------------------------------
# ----------------------------------------------------------------

# ------------------- 3D Compression Attenuation  ------------------
simulation_data = loadData3d("out/processed/3d_80Kby15_V3.jld2") # 15by15 tiles out to 300 long channel

plot_ωγ_attenuation_2d(simulation_data, .5, 1.2)    

# ------------------- 3D Sitched Attenuation  ------------------
gamma_values = 10 .^ range(-2, stop=2, length=5)
gamma_values = [ .1, .2, .5]
plotStitchAttenuation(simulation_data, gamma_values, 1.2)



# ----------------------------------------------------------------
# ------------------- 2D  -----------------------------------
# ----------------------------------------------------------------
# ------------------- 2D Compression MeanField  ------------------

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

# ------------------- 2D CompressionMode Density  ------------------
    mat"""plotDampedModeDensityPDF("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N400_K100_M1.mat", [.2, .01, .001], [.1])"""

# ------------------- 2D Shear Mean Field------------------

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