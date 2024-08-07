using MAT
using Glob
using Debugger

    
mutable struct ellipse_data
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
end

function main()
    simulation_data = ellipse_cruncher()
end

function plot_ellipse_pdf(simulation_data, gamma_value)
    closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
    matching_gamma_data = filter(entry -> entry.gamma == gamma_value, simulation_data)

    pressure_list = unique([entry.pressure for entry in matching_gamma_data]) # goes through each entry of simulation_data and get the P value at that entry
    
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    for pressure_value in pressure_list
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]
        
        # Each omega gamma value spans all seeds
        omega_gamma_list = unique([entry.omega_gamma for entry in matching_pressure_data])
        for omega_gamma_value in omega_gamma_list
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data)
            # mean(matrix, dims=2) means across rows. Dims = 1 is cross columns, dims =3 is into the screen
            #  If you have an array [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]], splatting converts this into calling hcat([a1, a2, a3], [b1, b2, b3], [c1, c2, c3]).
            mean_asp_rat_counts = mean(hcat([entry.asp_rat_counts for entry in matching_omega_gamma_data]...), dims=2)
            mean_asp_rat_bins = mean(hcat([entry.asp_rat_bins for entry in matching_omega_gamma_data]...), dims=2)
            mean_rot_ang_counts = mean(hcat([entry.rot_ang_counts for entry in matching_omega_gamma_data]...), dims=2)
            mean_rot_ang_bins = mean(hcat([entry.rot_ang_bins for entry in matching_omega_gamma_data]...), dims=2)

        end
        closest_asp_rat_counts = matching_pressure_data[closest_gamma_index].asp_rat_counts
    end


    closest_asp_rat_counts = simulation_data[closest_gamma_index].asp_rat_counts
end


function ellipse_cruncher()


    directory = "out/simulation_2d/bi_K100_all_seeds/"

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
            iloop_file_data["input_pressure"],
            iloop_file_data["driving_angular_frequency_dimensionless"],
            iloop_file_data["gamma_dimensionless"],
            asp_rat_counts,
            asp_rat_bins,
            rot_ang_counts,
            rot_ang_bins,
            iloop_file_data["driving_angular_frequency_dimensionless"]* iloop_file_data["gamma_dimensionless"],
            iloop_file_data["seed"],
            iloop_file_data["pressure_dimensionless"],
            -iloop_file_data["attenuation_x_dimensionless"],
            -iloop_file_data["attenuation_y_dimensionless"],
            iloop_file_data["wavespeed_x"],
            iloop_file_data["wavenumber_x"],
            iloop_file_data["mean_aspect_ratio"],
            iloop_file_data["mean_rotation_angles"],
            maximum(iloop_file_data["initial_distance_from_oscillation_output_x_fft"]),
            iloop_file_data["ellipse_stats_nonzero"],
            maximum(iloop_file_data["initial_distance_from_oscillation_output_y_fft"])
        )

        push!(simulation_data, data_entry)
    end

    return simulation_data
end