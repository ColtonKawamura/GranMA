# Example script for plotting eigenvalue data from MATLAB files

using MAT
using Plots
include("src/plotRealImagEigenValues.jl")

# Load the data file
vars = matread("out/2d_damped_eigenStuff/combined/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat")
outData = vars["outData"]

# Print available pressure and damping values to help with debugging
println("Available pressures: ", unique(round.(vec(outData["pressure"]), digits=5)))
println("Available damping values: ", unique(outData["damping"]))

# Basic plot example - low damping (using actual values from the dataset)
p1 = plotRealImagEigenValues(outData, [0.00118, 0.00574, 0.01388, 0.05882, 0.1, 0.2174], [0.001], loglog=true)
slopeLine("loglog", 0.25, [0.07, 0.4], 1.3, TextLocation=[0.1, 1.5])
slopeLine("loglog", 1, [0.5, 5], 2, TextLocation=[1, 3])
slopeLine("loglog", 2, [0.07, 0.4], 0.1, TextLocation=[0.08, 0.05])
savefig(p1, "julia_eigenvalues_low_damping_loglog.png")

# With scaling - low damping
p2 = plotRealImagEigenValues(outData, [0.00118, 0.00574, 0.01388, 0.05882, 0.1, 0.2174], [0.001], loglog=true, scaling=true)
slopeLine("loglog", 1, [2, 80], 15, TextLocation=[10, 20])
slopeLine("loglog", 0.25, [1, 9], 30, TextLocation=[3, 40])
savefig(p2, "julia_eigenvalues_low_damping_scaled_loglog.png")

# Medium damping (use a damping value that exists in the dataset)
p3 = plotRealImagEigenValues(outData, [0.00118, 0.00574, 0.01388, 0.05882, 0.1, 0.2174], [0.05], loglog=true)
slopeLine("loglog", 1, [0.5, 5], 2, TextLocation=[1, 3])
savefig(p3, "julia_eigenvalues_medium_damping_loglog.png")

# Medium damping with scaling
p4 = plotRealImagEigenValues(outData, [0.00118, 0.00574, 0.01388, 0.05882, 0.1, 0.2174], [0.05], loglog=true, scaling=true)
slopeLine("loglog", 1, [2, 80], 15, TextLocation=[10, 20])
slopeLine("loglog", 0.25, [1, 9], 30, TextLocation=[3, 40])
savefig(p4, "julia_eigenvalues_medium_damping_scaled_loglog.png")

# High damping (check if 0.25 is available, otherwise use a close value)
p5 = plotRealImagEigenValues(outData, [0.00118, 0.00574, 0.01388, 0.05882, 0.1, 0.2174], [0.25])
slopeLine("loglog", 2, [0.1, 2], 0.1, TextLocation=[0.5, 0.05])
savefig(p5, "julia_eigenvalues_high_damping_loglog.png")

println("All plots generated successfully!")


using GLMakie
using LaTeXStrings

f = Figure()
ax = Axis(f[1, 1],
    xlabel = L"Im(\lambda)/\sqrt{\hat{P}}",
    ylabel = L"Re(\lambda)/(\sqrt{\hat{P}\hat{\gamma}})",
    title = L"\text{My LaTeX Title}"
)

# Example with LaTeX tick labels
ax.xticks = ([0, 1, 2, 3], [L"0", L"1", L"2", L"3"])
ax.yticks = ([0, 0.5, 1.0], [L"0", L"\frac{1}{2}", L"1"])

scatter!(ax, rand(10), rand(10))
f


