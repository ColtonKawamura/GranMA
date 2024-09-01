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

export  makeJobList2D,
        load_data,
        file_data,
        crunch_and_save,
        crunch, 
        save_data,
        simulation_2d,
        pack_poly_2d, 
        plot_ellipse_low_pressure,
        process_outputs_2d,
        FilterData,
        para_crunch_and_save,
        para_crunch

mutable struct file_data
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
    wavenumber_y::Float64
    alphaoveromega_x::Float64
    alphaoveromega_y::Float64
    amplitude_vector_x::Vector{Float64}
    amplitude_vector_y::Vector{Float64}
    unwrapped_phase_vector_x::Vector{Float64}
    unwrapped_phase_vector_y::Vector{Float64}
    initial_distance_from_oscillation_output_x_fft::Vector{Float64}
    initial_distance_from_oscillation_output_y_fft::Vector{Float64}
    yfft_posY::Vector{Float64}
    yfft_posX::Vector{Float64}
    yfft_init_posX::Vector{Float64}
    yfft_init_posY::Vector{Float64}
    xfft_posY::Vector{Float64}
    xfft_posX::Vector{Float64}
    xfft_init_posX::Vector{Float64}
    xfft_init_pos::Vector{Float64}
    timeVec::Vector{Float64}
end

function makeJobList2D(filename::String, K_values::Vector{T1}, M_values::Vector{T2}, Bv_values::Vector{T3}, w_D_values::Vector{T4}, N_values::Vector{T5}, P_values::Vector{T6}, W_values::Vector{T7}, seeds::Vector{T8}) where {T1, T2, T3, T4, T5, T6, T7, T8}
    # makeJobList2D("testfile.txt", [100],[1],exp10.(-5:.2:2),exp10.(-5:.2:2),[10000],[.01,.001],[10,20,50],[1,2,3,4,5])
    # makeJobList2D("simulation_job_list.txt", [100],[1],exp10.(-5:.2:2),exp10.(-5:.2:2),[5000],[0.0001,0.00021544, 0.00046416, 0.001, 0.0021544, 0.0046416, 0.01, 0.021544, 0.046416, 0.1, ],[5],[1,2,3,4,5])    # Convert input vectors to the required types
    K_values = Int64.(K_values)
    M_values = Int64.(M_values)
    Bv_values = Float64.(Bv_values)
    w_D_values = Float64.(w_D_values)
    N_values = Float64.(N_values)
    P_values = Float64.(P_values)
    W_values = Int64.(W_values)
    seeds = Int64.(seeds)

    # Function to generate the MATLAB command EXAMPLE simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    function generate_matlab_command(K, M, Bv, w_D, N, P, W, seed)
        return "matlab -nodisplay -nosplash -r \"addpath('./src/matlab_functions/'); simulation_2d($K, $M, $Bv, $w_D, $N, $P, $W, $seed); exit\""
    end

    # Generate all combinations of parameters using IterTools.product
    combinations = IterTools.product(K_values, M_values, Bv_values, w_D_values, N_values, P_values, W_values, seeds)

    # Open a file to write the MATLAB commands
    open(filename, "w") do file
        for combo in combinations
            command = generate_matlab_command(combo...)
            println(file, command)
        end
    end

    println("Commands written to $filename")
end

function crunch_and_save(datapath::String, filepath::String)
    simulation_data = crunch(datapath)
    save_data(simulation_data, filepath)
    
    # Load the data back to verify
    reloaded_data = load_data(filepath)
    println("Data reloaded successfully. Number of entries: ", length(reloaded_data))
end

function crunch(datapath::String)
    # directory = "out/simulation_2d/bi_K100_all_seeds/"
    # directory = "out/simulation_2d/"
    mat_files = glob("*.mat", datapath)
    simulation_data = file_data[]
    
    for file_name in mat_files
        iloop_file_data = matread(file_name)

        # Extract data fields
        input_pressure = iloop_file_data["input_pressure"]
        omega = iloop_file_data["driving_angular_frequency_dimensionless"]
        gamma = iloop_file_data["gamma_dimensionless"]
        pressure_actual = iloop_file_data["pressure_dimensionless"]
        attenuation_x = iloop_file_data["attenuation_x_dimensionless"]
        attenuation_y = iloop_file_data["attenuation_y_dimensionless"]
        wavespeed_x = iloop_file_data["wavespeed_x"]
        wavenumber_x = iloop_file_data["wavenumber_x_dimensionless"]
        mean_aspect_ratio = iloop_file_data["mean_aspect_ratio"]
        mean_rotation_angles = iloop_file_data["mean_rotation_angles"]
        wavenumber_y = iloop_file_data["wavenumber_y_dimensionless"]

        # Check for NaN values
        if isnan(mean_rotation_angles) || isnan(mean_aspect_ratio) || 
           isnan(wavenumber_x) || isnan(wavenumber_y) ||
           isnan(wavespeed_x) || isnan(input_pressure) || 
           isnan(omega) || isnan(gamma) || isnan(pressure_actual) || 
           isnan(attenuation_x) || isnan(attenuation_y)
            println("Skipping file due to NaN values: $file_name")
            continue
        end

        # Convert matrix-type data to vectors
        asp_rat_counts = vec(iloop_file_data["asp_rat_counts"])
        asp_rat_bins = vec(iloop_file_data["asp_rat_bins"])
        rot_ang_counts = vec(iloop_file_data["rot_ang_counts"])
        rot_ang_bins = vec(iloop_file_data["rot_ang_bins"])
        # Extract amplitude and unwrapped phase vectors with checks
        amplitude_vector_x = isa(iloop_file_data["amplitude_vector_x"], Array) ? 
                             vec(iloop_file_data["amplitude_vector_x"]) : 
                             [iloop_file_data["amplitude_vector_x"]]

        amplitude_vector_y = isa(iloop_file_data["amplitude_vector_y"], Array) ? 
                             vec(iloop_file_data["amplitude_vector_y"]) : 
                             [iloop_file_data["amplitude_vector_y"]]

        unwrapped_phase_vector_x = haskey(iloop_file_data, "unwrapped_phase_vector_x") ? 
                                   (isa(iloop_file_data["unwrapped_phase_vector_x"], Array) ? 
                                       vec(iloop_file_data["unwrapped_phase_vector_x"]) : 
                                       [iloop_file_data["unwrapped_phase_vector_x"]]) : 
                                   Float64[]

        unwrapped_phase_vector_y = haskey(iloop_file_data, "unwrapped_phase_vector_y") ? 
                                   (isa(iloop_file_data["unwrapped_phase_vector_y"], Array) ? 
                                       vec(iloop_file_data["unwrapped_phase_vector_y"]) : 
                                       [iloop_file_data["unwrapped_phase_vector_y"]]) : 
                                   Float64[]

        initial_distance_from_oscillation_output_x_fft = haskey(iloop_file_data, "initial_distance_from_oscillation_output_x_fft") ? 
                                                        (isa(iloop_file_data["initial_distance_from_oscillation_output_x_fft"], Array) ? 
                                                            vec(iloop_file_data["initial_distance_from_oscillation_output_x_fft"]) : 
                                                            [iloop_file_data["initial_distance_from_oscillation_output_x_fft"]]) : 
                                                        Float64[]

        initial_distance_from_oscillation_output_y_fft = haskey(iloop_file_data, "initial_distance_from_oscillation_output_y_fft") ? 
                                                        (isa(iloop_file_data["initial_distance_from_oscillation_output_y_fft"], Array) ? 
                                                            vec(iloop_file_data["initial_distance_from_oscillation_output_y_fft"]) : 
                                                            [iloop_file_data["initial_distance_from_oscillation_output_y_fft"]]) : 
                                                        Float64[]

        
        yfft_posY = isa(iloop_file_data["y_all_yfft"], Array) ? 
                            vec(iloop_file_data["y_all_yfft"]) : 
                            [iloop_file_data["y_all_yfft"]]

        yfft_posX = isa(iloop_file_data["x_all_yfft"], Array) ? 
                            vec(iloop_file_data["x_all_yfft"]) : 
                            [iloop_file_data["x_all_yfft"]]

        yfft_init_posX = isa(iloop_file_data["x0_yfft"], Array) ? 
                            vec(iloop_file_data["x0_yfft"]) : 
                            [iloop_file_data["x0_yfft"]]

        yfft_init_posY = isa(iloop_file_data["y0_yfft"], Array) ? 
                            vec(iloop_file_data["y0_yfft"]) : 
                            [iloop_file_data["y0_yfft"]]

        xfft_posY = isa(iloop_file_data["y_all_xfft"], Array) ? 
                            vec(iloop_file_data["y_all_xfft"]) : 
                            [iloop_file_data["y_all_xfft"]]

        xfft_posX = isa(iloop_file_data["x_all_xfft"], Array) ? 
                            vec(iloop_file_data["x_all_xfft"]) : 
                            [iloop_file_data["x_all_xfft"]]

        xfft_init_posX = isa(iloop_file_data["x0_xfft"], Array) ? 
                            vec(iloop_file_data["x0_xfft"]) : 
                            [iloop_file_data["x0_xfft"]]

        xfft_init_posY = isa(iloop_file_data["y0_xfft"], Array) ? 
                            vec(iloop_file_data["y0_xfft"]) : 
                            [iloop_file_data["y0_xfft"]]
        
        timeVec = isa(iloop_file_data["tvec"], Array) ? 
                            vec(iloop_file_data["tvec"]) : 
                            [iloop_file_data["tvec"]]


        # Handle potentially empty arrays for fft limits
        fft_x = get(iloop_file_data, "initial_distance_from_oscillation_output_x_fft", [])
        fft_y = get(iloop_file_data, "initial_distance_from_oscillation_output_y_fft", [])
        fft_limit_x = isempty(fft_x) ? NaN : maximum(fft_x)
        fft_limit_y = isempty(fft_y) ? NaN : maximum(fft_y)

        data_entry = file_data(
            input_pressure,
            omega,
            gamma,
            asp_rat_counts,
            asp_rat_bins,
            rot_ang_counts,
            rot_ang_bins,
            omega * gamma,
            iloop_file_data["seed"],
            pressure_actual,
            -attenuation_x,
            -attenuation_y,
            -wavespeed_x,
            -wavenumber_x,
            mean_aspect_ratio,
            mean_rotation_angles,
            fft_limit_x,
            iloop_file_data["ellipse_stats_nonzero"],
            fft_limit_y,
            wavenumber_y,
            -attenuation_x / omega,
            -attenuation_y / omega,
            amplitude_vector_x,
            amplitude_vector_y,
            unwrapped_phase_vector_x,
            unwrapped_phase_vector_y,
            initial_distance_from_oscillation_output_x_fft,
            initial_distance_from_oscillation_output_y_fft,
            yfft_posY,
            yfft_posX,
            yfft_init_posX,
            yfft_init_posY,
            xfft_posY,
            xfft_posX,
            xfft_init_posX,
            xfft_init_posY,
            timeVec
        )

        push!(simulation_data, data_entry)
    end

    return simulation_data
end

function save_data(simulation_data::Vector{file_data}, filepath::String)
    # filepath = "out/processed/name_without_extension
    @save filepath simulation_data
    println("Data saved to $filepath")
end

function load_data(filepath::String)::Vector{file_data}
    # filename = "out/processed/name_without_extension
    filepath
    @load filepath simulation_data
    return simulation_data
end

function pack_poly_2d(N, K, D, G, M, P_target, W_factor, seed, plotit)
    mat"""
    addpath('src/matlab_functions/')
    N =double($(N))
    K = double($(K))
    M = double($(M))
    D = double($(D))
    G = double($(G))
    P_target = double($(P_target))
    W_factor = double($(W_factor))
    seed = double($(seed))
    plotit = double($(plotit))
    packing_poly_2d(N, K, D, G, M, P_target, W_factor, seed, plotit)
    """
end

function simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    mat"""
    addpath('src/matlab_functions/')
    K = double($(K))
    M = double($(M))
    Bv = double($(Bv))
    w_D = double($(w_D))
    N =double($(N))
    P = double($(P))
    W = double($(W))
    seed = double($(seed))
    simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    """
end

function FilterData(data::Vector{file_data}, args...)
    # Initialize the filtered data to be the full simulation_data
    filtered_data = data

    # Iterate over the provided arguments in pairs
    for i in 1:2:length(args)
        value = args[i]
        field_name = args[i+1]

        # Determine the closest match for the given field
        closest_index = argmin(abs.([getfield(idx, field_name) for idx in filtered_data] .- value))
        closest_value = getfield(filtered_data[closest_index], field_name)

        # Filter the data to match the closest value
        filtered_data = filter(entry -> getfield(entry, field_name) == closest_value, filtered_data)
    end

    return filtered_data
end

function para_crunch(datapath::String)
    # directory = "out/simulation_2d/bi_K100_all_seeds/"
    # directory = "out/simulation_2d/"
    mat_files = glob("*.mat", datapath)
    thread_data = [Vector{file_data}() for _ in 1:Threads.nthreads()]
    
    Threads.@threads for file_name in mat_files
        iloop_file_data = matread(file_name)

        # Extract data fields
        input_pressure = iloop_file_data["input_pressure"]
        omega = iloop_file_data["driving_angular_frequency_dimensionless"]
        gamma = iloop_file_data["gamma_dimensionless"]
        pressure_actual = iloop_file_data["pressure_dimensionless"]
        attenuation_x = iloop_file_data["attenuation_x_dimensionless"]
        attenuation_y = iloop_file_data["attenuation_y_dimensionless"]
        wavespeed_x = iloop_file_data["wavespeed_x"]
        wavenumber_x = iloop_file_data["wavenumber_x_dimensionless"]
        mean_aspect_ratio = iloop_file_data["mean_aspect_ratio"]
        mean_rotation_angles = iloop_file_data["mean_rotation_angles"]
        wavenumber_y = iloop_file_data["wavenumber_y_dimensionless"]

        # Check for NaN values
        if isnan(mean_rotation_angles) || isnan(mean_aspect_ratio) || 
           isnan(wavenumber_x) || isnan(wavenumber_y) ||
           isnan(wavespeed_x) || isnan(input_pressure) || 
           isnan(omega) || isnan(gamma) || isnan(pressure_actual) || 
           isnan(attenuation_x) || isnan(attenuation_y)
            println("Skipping file due to NaN values: $file_name")
            continue
        end

        # Convert matrix-type data to vectors
        asp_rat_counts = vec(iloop_file_data["asp_rat_counts"])
        asp_rat_bins = vec(iloop_file_data["asp_rat_bins"])
        rot_ang_counts = vec(iloop_file_data["rot_ang_counts"])
        rot_ang_bins = vec(iloop_file_data["rot_ang_bins"])
        # Extract amplitude and unwrapped phase vectors with checks
        amplitude_vector_x = isa(iloop_file_data["amplitude_vector_x"], Array) ? 
                             vec(iloop_file_data["amplitude_vector_x"]) : 
                             [iloop_file_data["amplitude_vector_x"]]

        amplitude_vector_y = isa(iloop_file_data["amplitude_vector_y"], Array) ? 
                             vec(iloop_file_data["amplitude_vector_y"]) : 
                             [iloop_file_data["amplitude_vector_y"]]

        unwrapped_phase_vector_x = haskey(iloop_file_data, "unwrapped_phase_vector_x") ? 
                                   (isa(iloop_file_data["unwrapped_phase_vector_x"], Array) ? 
                                       vec(iloop_file_data["unwrapped_phase_vector_x"]) : 
                                       [iloop_file_data["unwrapped_phase_vector_x"]]) : 
                                   Float64[]

        unwrapped_phase_vector_y = haskey(iloop_file_data, "unwrapped_phase_vector_y") ? 
                                   (isa(iloop_file_data["unwrapped_phase_vector_y"], Array) ? 
                                       vec(iloop_file_data["unwrapped_phase_vector_y"]) : 
                                       [iloop_file_data["unwrapped_phase_vector_y"]]) : 
                                   Float64[]

        initial_distance_from_oscillation_output_x_fft = haskey(iloop_file_data, "initial_distance_from_oscillation_output_x_fft") ? 
                                                        (isa(iloop_file_data["initial_distance_from_oscillation_output_x_fft"], Array) ? 
                                                            vec(iloop_file_data["initial_distance_from_oscillation_output_x_fft"]) : 
                                                            [iloop_file_data["initial_distance_from_oscillation_output_x_fft"]]) : 
                                                        Float64[]

        initial_distance_from_oscillation_output_y_fft = haskey(iloop_file_data, "initial_distance_from_oscillation_output_y_fft") ? 
                                                        (isa(iloop_file_data["initial_distance_from_oscillation_output_y_fft"], Array) ? 
                                                            vec(iloop_file_data["initial_distance_from_oscillation_output_y_fft"]) : 
                                                            [iloop_file_data["initial_distance_from_oscillation_output_y_fft"]]) : 
                                                        Float64[]


        # Handle potentially empty arrays for fft limits
        fft_x = get(iloop_file_data, "initial_distance_from_oscillation_output_x_fft", [])
        fft_y = get(iloop_file_data, "initial_distance_from_oscillation_output_y_fft", [])
        fft_limit_x = isempty(fft_x) ? NaN : maximum(fft_x)
        fft_limit_y = isempty(fft_y) ? NaN : maximum(fft_y)

        data_entry = file_data(
            input_pressure,
            omega,
            gamma,
            asp_rat_counts,
            asp_rat_bins,
            rot_ang_counts,
            rot_ang_bins,
            omega * gamma,
            iloop_file_data["seed"],
            pressure_actual,
            -attenuation_x,
            -attenuation_y,
            -wavespeed_x,
            -wavenumber_x,
            mean_aspect_ratio,
            mean_rotation_angles,
            fft_limit_x,
            iloop_file_data["ellipse_stats_nonzero"],
            fft_limit_y,
            wavenumber_y,
            -attenuation_x / omega,
            -attenuation_y / omega,
            amplitude_vector_x,
            amplitude_vector_y,
            unwrapped_phase_vector_x,
            unwrapped_phase_vector_y,
            initial_distance_from_oscillation_output_x_fft,
            initial_distance_from_oscillation_output_y_fft
        )

        push!(thread_data[Threads.threadid()], data_entry) # pushes to respective threads
    end

    # Combine results from all threads
    simulation_data = vcat(thread_data...)
    return simulation_data
end

function para_crunch_and_save(datapath::String, filepath::String)
    simulation_data = para_crunch(datapath)
    save_data(simulation_data, filepath)
    
    # Load the data back to verify
    reloaded_data = load_data(filepath)
    println("Data reloaded successfully. Number of entries: ", length(reloaded_data))
end

end