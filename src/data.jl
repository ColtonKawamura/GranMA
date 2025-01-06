include("types.jl")
export 
    crunch_and_save,
    crunch,
    save_data,
    load_data,
    FilterData,
    para_crunch,
    para_crunch_and_save,
    paraSimpleCrunch,
    crunchGausData,
    saveGausData,
    loadGausData,
    crunchNSaveGaus,
    filterDataGaus,
    crunch3d

function filterDataGaus(data::Vector{gaus_data}, args...)
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

function crunchNSaveGaus(datapath::String, filepath::String)
    simulation_data = crunchGausData(datapath)
    saveGausData(simulation_data, filepath)
    
    # Load the data back to verify
    reloaded_data = loadGausData(filepath)
    println("Data reloaded successfully. Number of entries: ", length(reloaded_data))
end

function loadGausData(filepath::String)::Vector{gaus_data}
    # filename = "out/processed/name_without_extension
    filepath
    @load filepath simulation_data
    return simulation_data
end

function saveGausData(simulation_data::Vector{gaus_data}, filepath::String)
     # filepath = "out/processed/name_without_extension
    @save filepath simulation_data
    println("Data saved to $filepath")
end
function crunchGausData(datapath::String)
    mat_files = glob("*.mat", datapath)
    simulation_data = gaus_data[]
    for file_name in mat_files
        iloop_file_data = matread(file_name)
        
        # Extract data from the .mat file
        mean_diameter = iloop_file_data["meanDiameter"]
        omega = iloop_file_data["omega"]
        gamma = iloop_file_data["gamma"]
        spring_constant = iloop_file_data["spring_constant"]
        mass = iloop_file_data["mass"]
        pressure = iloop_file_data["pressure"]
        width = iloop_file_data["width"]
        seed = iloop_file_data["seed"]
        pressure_actual = iloop_file_data["pressure_actual"]
        attenuation = iloop_file_data["attenuation"]
        wavespeed = iloop_file_data["wavespeed"]
        wavenumber = iloop_file_data["wavenumber"]
        dt = iloop_file_data["dt"]

        # Add data to the simulation_data vector as gaus_data struct
        push!(simulation_data, gaus_data(mean_diameter, omega, gamma, spring_constant, mass, pressure, width, seed, pressure_actual, attenuation, wavespeed, wavenumber, dt))
    end
    return simulation_data
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

        
        # yfft_posY = isa(iloop_file_data["y_all_yfft"], Array) ? 
        #                     vec(iloop_file_data["y_all_yfft"]) : 
        #                     [iloop_file_data["y_all_yfft"]]

        # yfft_posX = isa(iloop_file_data["x_all_yfft"], Array) ? 
        #                     vec(iloop_file_data["x_all_yfft"]) : 
        #                     [iloop_file_data["x_all_yfft"]]

        # yfft_init_posX = isa(iloop_file_data["x0_yfft"], Array) ? 
        #                     vec(iloop_file_data["x0_yfft"]) : 
        #                     [iloop_file_data["x0_yfft"]]

        # yfft_init_posY = isa(iloop_file_data["y0_yfft"], Array) ? 
        #                     vec(iloop_file_data["y0_yfft"]) : 
        #                     [iloop_file_data["y0_yfft"]]

        # xfft_posY = isa(iloop_file_data["y_all_xfft"], Array) ? 
        #                     vec(iloop_file_data["y_all_xfft"]) : 
        #                     [iloop_file_data["y_all_xfft"]]

        # xfft_posX = isa(iloop_file_data["x_all_xfft"], Array) ? 
        #                     vec(iloop_file_data["x_all_xfft"]) : 
        #                     [iloop_file_data["x_all_xfft"]]

        # xfft_init_posX = isa(iloop_file_data["x0_xfft"], Array) ? 
        #                     vec(iloop_file_data["x0_xfft"]) : 
        #                     [iloop_file_data["x0_xfft"]]

        # xfft_init_posY = isa(iloop_file_data["y0_xfft"], Array) ? 
        #                     vec(iloop_file_data["y0_xfft"]) : 
        #                     [iloop_file_data["y0_xfft"]]
        
        # timeVec = isa(iloop_file_data["tvec"], Array) ? 
        #                     vec(iloop_file_data["tvec"]) : 
        #                     [iloop_file_data["tvec"]]


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
            # yfft_posY,
            # yfft_posX,
            # yfft_init_posX,
            # yfft_init_posY,
            # xfft_posY,
            # xfft_posX,
            # xfft_init_posX,
            # xfft_init_posY,
            # timeVec
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


        # yfft_posY = isa(iloop_file_data["y_all_yfft"], Array) ? 
        #                     vec(iloop_file_data["y_all_yfft"]) : 
        #                     [iloop_file_data["y_all_yfft"]]
        
        # yfft_posX = isa(iloop_file_data["x_all_yfft"], Array) ? 
        #                     vec(iloop_file_data["x_all_yfft"]) : 
        #                     [iloop_file_data["x_all_yfft"]]

        # yfft_init_posX = isa(iloop_file_data["x0_yfft"], Array) ? 
        #                     vec(iloop_file_data["x0_yfft"]) : 
        #                     [iloop_file_data["x0_yfft"]]

        # yfft_init_posY = isa(iloop_file_data["y0_yfft"], Array) ? 
        #                     vec(iloop_file_data["y0_yfft"]) : 
        #                     [iloop_file_data["y0_yfft"]]

        # xfft_posY = isa(iloop_file_data["y_all_xfft"], Array) ? 
        #                     vec(iloop_file_data["y_all_xfft"]) : 
        #                     [iloop_file_data["y_all_xfft"]]

        # xfft_posX = isa(iloop_file_data["x_all_xfft"], Array) ? 
        #                     vec(iloop_file_data["x_all_xfft"]) : 
        #                     [iloop_file_data["x_all_xfft"]]

        # xfft_init_posX = isa(iloop_file_data["x0_xfft"], Array) ? 
        #                     vec(iloop_file_data["x0_xfft"]) : 
        #                     [iloop_file_data["x0_xfft"]]

        # xfft_init_posY = isa(iloop_file_data["y0_xfft"], Array) ? 
        #                     vec(iloop_file_data["y0_xfft"]) : 
        #                     [iloop_file_data["y0_xfft"]]
        
        # timeVec = isa(iloop_file_data["tvec"], Array) ? 
        #                     vec(iloop_file_data["tvec"]) : 
        #                     [iloop_file_data["tvec"]]


        # # Handle potentially empty arrays for fft limits
        # fft_x = get(iloop_file_data, "initial_distance_from_oscillation_output_x_fft", [])
        # fft_y = get(iloop_file_data, "initial_distance_from_oscillation_output_y_fft", [])
        # fft_limit_x = isempty(fft_x) ? NaN : maximum(fft_x)
        # fft_limit_y = isempty(fft_y) ? NaN : maximum(fft_y)

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
            # yfft_posY,
            # yfft_posX,
            # yfft_init_posX,
            # yfft_init_posY,
            # xfft_posY,
            # xfft_posX,
            # xfft_init_posX,
            # xfft_init_posY,
            # timeVec
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

function paraSimpleCrunch(datapath::String)
    mat_files = glob("*.mat", datapath)
    thread_data = [Vector{file_data}() for _ in 1:Threads.nthreads()]
    
    Threads.@threads for file_name in mat_files
        iloop_file_data = matread(file_name)

        # Extract data fields
        width = iloop_file_data["W"]
        pressureActual = iloop_file_data["P"]
        dt = iloop_file_data["dt"]
        particles = size(iloop_file_data["Pos2D"], 2)  # Extract the number of particles
        timesteps = Int(iloop_file_data["Nt"])
        mass = iloop_file_data["M"]
        amplitude = iloop_file_data["A"]
        pressureInput = iloop_file_data["input_pressure"]
        gamma = iloop_file_data["Bv"]
        omega = iloop_file_data["w_D"]
        seed = iloop_file_data["seed"]
        kappa = iloop_file_data["K"]
        diameters = iloop_file_data["diameters"]

        trajectory = [
            ([Pos2D(iloop_file_data["Pos2D"][t, p, 1], iloop_file_data["Pos2D"][t, p, 2]) for p in 1:particles], diameters)
            for t in 1:timesteps
        ]

        data_entry = raw_data(
            width,
            pressureActual,
            dt,
            trajectory,
            timesteps,
            mass,
            amplitude,
            pressureInput,
            gamma/sqrt(kappa*mass),
            omega*sqrt(mass/kappa),
            seed,
            kappa,
            gamma/sqrt(kappa*mass) * omega*sqrt(mass/kappa),
        )

        push!(thread_data[Threads.threadid()], data_entry) # pushes to respective threads
    end

    # Combine results from all threads
    simulation_data = vcat(thread_data...)
    return simulation_data
end


function crunch3d(datapath::String)
    # directory = "out/simulation_3d/initial_test/"
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
        attenuation_z = iloop_file_data["attenuation_z_dimensionless"]
        wavespeed_x = iloop_file_data["wavespeed_x"]
        wavenumber_x = iloop_file_data["wavenumber_x_dimensionless"]
        wavenumber_y = iloop_file_data["wavenumber_y_dimensionless"]
        wavenumber_z = iloop_file_data["wavenumber_z_dimensionless"]

        # Check for NaN values
        if isnan(mean_rotation_angles) || isnan(mean_aspect_ratio) || 
           isnan(wavenumber_x) || isnan(wavenumber_y) ||
           isnan(wavespeed_x) || isnan(input_pressure) || 
           isnan(omega) || isnan(gamma) || isnan(pressure_actual) || 
           isnan(attenuation_x) || isnan(attenuation_y)
            println("Skipping file due to NaN values: $file_name")
            continue
        end

        # Extract amplitude and unwrapped phase vectors with checks
        amplitude_vector_x = isa(iloop_file_data["amplitude_vector_x"], Array) ? 
                             vec(iloop_file_data["amplitude_vector_x"]) : 
                             [iloop_file_data["amplitude_vector_x"]]

        amplitude_vector_y = isa(iloop_file_data["amplitude_vector_y"], Array) ? 
                             vec(iloop_file_data["amplitude_vector_y"]) : 
                             [iloop_file_data["amplitude_vector_y"]]

        amplitude_vector_z = isa(iloop_file_data["amplitude_vector_z"], Array) ? 
                             vec(iloop_file_data["amplitude_vector_z"]) : 
                             [iloop_file_data["amplitude_vector_z"]]

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

        unwrapped_phase_vector_z = haskey(iloop_file_data, "unwrapped_phase_vector_z") ? 
                                   (isa(iloop_file_data["unwrapped_phase_vector_z"], Array) ? 
                                       vec(iloop_file_data["unwrapped_phase_vector_z"]) : 
                                       [iloop_file_data["unwrapped_phase_vector_z"]]) : 
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

        initial_distance_from_oscillation_output_z_fft = haskey(iloop_file_data, "initial_distance_from_oscillation_output_z_fft") ? 
                                                        (isa(iloop_file_data["initial_distance_from_oscillation_output_z_fft"], Array) ? 
                                                            vec(iloop_file_data["initial_distance_from_oscillation_output_z_fft"]) : 
                                                            [iloop_file_data["initial_distance_from_oscillation_output_z_fft"]]) : 
                                                        Float64[]
        

        # Handle potentially empty arrays for fft limits
        fft_x = get(iloop_file_data, "initial_distance_from_oscillation_output_x_fft", [])
        fft_y = get(iloop_file_data, "initial_distance_from_oscillation_output_y_fft", [])
        fft_z = get(iloop_file_data, "initial_distance_from_oscillation_output_z_fft", [])
        fft_limit_x = isempty(fft_x) ? NaN : maximum(fft_x)
        fft_limit_y = isempty(fft_y) ? NaN : maximum(fft_y)
        fft_limit_z = isempty(fft_z) ? NaN : maximum(fft_z)

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
            initial_distance_from_oscillation_output_z_fft
       )

        push!(simulation_data, data_entry)
    end

    return simulation_data
end