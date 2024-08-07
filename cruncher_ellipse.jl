using MAT
using Glob
using Debugger

    
mutable struct ellipse_data
    P::Float64
    omega::Float64
    gamma::Float64
    asp_rat_counts::Vector{Float64}
    asp_rat_bins::Vector{Float64}
    rot_ang_counts::Vector{Float64}
    rot_ang_bins::Vector{Float64}
end

function main()
    simulation_data = ellipse_cruncher()
end

function plot_ellipse_pdf(simulation_data, gamma_value)
    closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
    matching_gamma_data = filter(entry -> entry.gamma == gamma_value, simulation_data)

    pressure_list = unique([entry.P for entry in matching_gamma_data]) # goes through each entry of simulation_data and get the P value at that entry

    for pressure_value in pressure_list
        matching_pressure_data = filter(entry -> entry.P == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
        closest_asp_rat_counts = matching_pressure_data[closest_gamma_index].asp_rat_counts
    end


    closest_asp_rat_counts = simulation_data[closest_gamma_index].asp_rat_counts
end


function ellipse_cruncher()


    directory = "out/simulation_2d/"

    mat_files = glob("*.mat", directory)

    simulation_data = ellipse_data[]

    for file_name in mat_files

        iloop_file_data = matread(file_name)
        
        # Need to do this because MATALB stored these as a matrix-type
        asp_rat_counts = vec(iloop_file_data["asp_rat_counts"])
        asp_rat_bins = vec(iloop_file_data["asp_rat_bins"])
        rot_ang_counts = vec(iloop_file_data["rot_ang_counts"])
        rot_ang_bins = vec(iloop_file_data["rot_ang_bins"])
    

        data_entry = ellipse_data(
            iloop_file_data["pressure_dimensionless"],
            iloop_file_data["driving_angular_frequency_dimensionless"],
            iloop_file_data["gamma_dimensionless"],
            asp_rat_counts,
            asp_rat_bins,
            rot_ang_counts,
            rot_ang_bins
        )

        push!(simulation_data, data_entry)
    end

    return simulation_data
end