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
using Distributed
using StaticArrays
using DSP
using Polynomials

include("jobs.jl")
include("data.jl")
include("simulation.jl")
include("types.jl")
include("plotting.jl")

end