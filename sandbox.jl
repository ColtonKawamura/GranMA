include("src/GranMA.jl")

using .GranMA
using CSV
using DataFrames

# Read the data table
data_frame = CSV.read("out/processed_simulation/K100_ellipse_edits.csv", DataFrame)

plot_ellipse_ωγ_2d(data_frame, .5,1)

simulation_2d(100, 1, .1, 1, 5000, .1, 5, 5)

plot_ωγ_attenuation_2d(data_frame, .1, 1)

plot_ωγ_wavespeed_2d(data_frame, .1)