include("src/GranMA.jl")

using .GranMA
using CSV
using DataFrames

# Read the data table
data_frame = CSV.read("out/test_crunch.csv", DataFrame)

plot_ellipse_pressure_2d(data_frame, .1)

simulation_2d(100, 1, .1, 1, 5000, .1, 5, 5)