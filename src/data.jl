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
    crunch3d,
    loadData3d,
    saveData3d,
    crunchNSave3d,
    FilterData3d

function FilterData3d(data::Vector{data3d}, args...)

    # Iterate over the provided arguments in pairs
    for i in 1:2:length(args)
        if isa(args[i], Vector)
            limits = args[i]
            dimension = args[i+1]
            if dimension == "y"
                lower_y = limits[1]
                upper_y = limits[2]
                for j in 1:length(data)
                    yFFT_initialY = data[j].y_fft_initial_y
                    mask = (yFFT_initialY .>= lower_y) .& (yFFT_initialY .<= upper_y)
                    data[j].y_fft_initial_y = data[j].y_fft_initial_y[mask]
                    data[j].y_fft_initial_z = data[j].y_fft_initial_z[mask]
                    data[j].initial_distance_from_oscillation_output_y_fft = data[j].initial_distance_from_oscillation_output_y_fft[mask]
                    data[j].amplitude_vector_y = data[j].amplitude_vector_y[mask]
                    data[j].unwrapped_phase_vector_y = data[j].unwrapped_phase_vector_y[mask]
                    
                    # This section filters out particles that have a initial z outside of [lower, upper] 
                    zFFT_initialY = data[j].z_fft_initial_y
                    mask = (zFFT_initialY .>= lower_y) .& (zFFT_initialY .<= upper_y)
                    data[j].z_fft_initial_y = data[j].z_fft_initial_y[mask]
                    data[j].z_fft_initial_z = data[j].z_fft_initial_z[mask]
                    data[j].initial_distance_from_oscillation_output_z_fft = data[j].initial_distance_from_oscillation_output_z_fft[mask]
                    data[j].amplitude_vector_z = data[j].amplitude_vector_z[mask]
                    data[j].unwrapped_phase_vector_z = data[j].unwrapped_phase_vector_z[mask]
    
                    # This section filters out particles that have a initial z outside of [lower, upper] for the x fft
                    xFFT_initialY = data[j].x_fft_initial_y
                    mask = (xFFT_initialY .>= lower_y) .& (xFFT_initialY .<= upper_y)
                    data[j].x_fft_initial_y = data[j].x_fft_initial_y[mask]
                    data[j].x_fft_initial_z = data[j].x_fft_initial_z[mask]
                    data[j].initial_distance_from_oscillation_output_x_fft = data[j].initial_distance_from_oscillation_output_x_fft[mask]
                    data[j].amplitude_vector_x = data[j].amplitude_vector_x[mask]
                    data[j].unwrapped_phase_vector_x = data[j].unwrapped_phase_vector_x[mask]
                end
            elseif dimension == "z"
                lower_z = limits[1]
                upper_z = limits[2]
                for j in 1:length(data)
                    yFFT_initialZ = data[j].y_fft_initial_z
                    mask = (yFFT_initialZ .>= lower_z) .& (yFFT_initialZ .<= upper_z)
                    data[j].y_fft_initial_y = data[j].y_fft_initial_y[mask]
                    data[j].y_fft_initial_z = data[j].y_fft_initial_z[mask]
                    data[j].initial_distance_from_oscillation_output_y_fft = data[j].initial_distance_from_oscillation_output_y_fft[mask]
                    data[j].amplitude_vector_y = data[j].amplitude_vector_y[mask]
                    data[j].unwrapped_phase_vector_y = data[j].unwrapped_phase_vector_y[mask]
        
                    zFFT_initialZ = data[j].z_fft_initial_z
                    mask = (zFFT_initialZ .>= lower_z) .& (zFFT_initialZ .<= upper_z)
                    data[j].z_fft_initial_y = data[j].z_fft_initial_y[mask]
                    data[j].z_fft_initial_z = data[j].z_fft_initial_z[mask]
                    data[j].initial_distance_from_oscillation_output_z_fft = data[j].initial_distance_from_oscillation_output_z_fft[mask]
                    data[j].amplitude_vector_z = data[j].amplitude_vector_z[mask]
                    data[j].unwrapped_phase_vector_z = data[j].unwrapped_phase_vector_z[mask]
        
                    xFFT_initialZ = data[j].x_fft_initial_z
                    mask = (xFFT_initialZ .>= lower_z) .& (xFFT_initialZ .<= upper_z)
                    data[j].x_fft_initial_y = data[j].x_fft_initial_y[mask]
                    data[j].x_fft_initial_z = data[j].x_fft_initial_z[mask]
                    data[j].initial_distance_from_oscillation_output_x_fft = data[j].initial_distance_from_oscillation_output_x_fft[mask]
                    data[j].amplitude_vector_x = data[j].amplitude_vector_x[mask]
                    data[j].unwrapped_phase_vector_x = data[j].unwrapped_phase_vector_x[mask]
                end
            else
                error("Invalid dimension specified. Use 'y' or 'z'.")
            end


        else
            value = args[i]
            field_name = args[i+1]
            # Determine the closest match for the given field
            closest_index = argmin(abs.([getfield(idx, field_name) for idx in data] .- value))
            closest_value = getfield(data[closest_index], field_name)

            # Filter the data to match the closest value
            data = filter(entry -> getfield(entry, field_name) == closest_value, data)
        end
    end

    return data
end

function saveData3d(simulation_data::Vector{data3d}, filepath::String)
    # filepath = "out/processed/name_without_extension
    @save filepath simulation_data
    println("Data saved to $filepath")
end

function loadData3d(filepath::String)::Vector{data3d}
    # filename = "out/processed/name_without_extension
    filepath
    @load filepath simulation_data
    return simulation_data
end

function crunchNSave3d(datapath::String, filepath::String)
    simulation_data = crunch3d(datapath)
    saveData3d(simulation_data, filepath)
    
    # Load the data back to verify
    reloaded_data = loadData3d(filepath)
    println("Data reloaded successfully. Number of entries: ", length(reloaded_data))
end

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
    simulation_data = data3d[]
    
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
        if isnan(wavenumber_x) || isnan(wavenumber_y) ||
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

        x_fft_initial_y = haskey(iloop_file_data, "x_fft_initial_y") ? 
                    (isa(iloop_file_data["x_fft_initial_y"], Array) ? 
                        vec(iloop_file_data["x_fft_initial_y"]) : 
                        [iloop_file_data["x_fft_initial_y"]]) : 
                    Float64[]

        x_fft_initial_z = haskey(iloop_file_data, "x_fft_initial_z") ?
                    (isa(iloop_file_data["x_fft_initial_z"], Array) ? 
                        vec(iloop_file_data["x_fft_initial_z"]) : 
                        [iloop_file_data["x_fft_initial_z"]]) : 
                    Float64[]
        
        y_fft_initial_y = haskey(iloop_file_data, "y_fft_initial_y") ?
                    (isa(iloop_file_data["y_fft_initial_y"], Array) ? 
                        vec(iloop_file_data["y_fft_initial_y"]) : 
                        [iloop_file_data["y_fft_initial_y"]]) : 
                    Float64[]

        y_fft_initial_z = haskey(iloop_file_data, "y_fft_initial_z") ?
                    (isa(iloop_file_data["y_fft_initial_z"], Array) ? 
                        vec(iloop_file_data["y_fft_initial_z"]) : 
                        [iloop_file_data["y_fft_initial_z"]]) : 
                    Float64[]
        
        z_fft_initial_y = haskey(iloop_file_data, "z_fft_initial_y") ?
                    (isa(iloop_file_data["z_fft_initial_y"], Array) ? 
                        vec(iloop_file_data["z_fft_initial_y"]) : 
                        [iloop_file_data["z_fft_initial_y"]]) : 
                    Float64[]
            
        z_fft_initial_z = haskey(iloop_file_data, "z_fft_initial_z") ?
                    (isa(iloop_file_data["z_fft_initial_z"], Array) ? 
                        vec(iloop_file_data["z_fft_initial_z"]) : 
                        [iloop_file_data["z_fft_initial_z"]]) : 
                    Float64[]
        

        # Handle potentially empty arrays for fft limits
        fft_x = get(iloop_file_data, "initial_distance_from_oscillation_output_x_fft", [])
        fft_y = get(iloop_file_data, "initial_distance_from_oscillation_output_y_fft", [])
        fft_z = get(iloop_file_data, "initial_distance_from_oscillation_output_z_fft", [])
        fft_limit_x = isempty(fft_x) ? NaN : maximum(fft_x)
        fft_limit_y = isempty(fft_y) ? NaN : maximum(fft_y)
        fft_limit_z = isempty(fft_z) ? NaN : maximum(fft_z)

        data_entry = data3d(
            input_pressure,
            omega,
            gamma,
            omega * gamma,
            iloop_file_data["seed"],
            pressure_actual,
            -attenuation_x,
            -attenuation_y,
            -attenuation_z,
            -wavespeed_x,
            -wavenumber_x,
            fft_limit_x,
            fft_limit_y,
            fft_limit_z,
            wavenumber_y,
            wavenumber_z,
            -attenuation_x / omega,
            -attenuation_y / omega,
            -attenuation_z / omega,
            amplitude_vector_x,
            amplitude_vector_y,
            amplitude_vector_z,
            unwrapped_phase_vector_x,
            unwrapped_phase_vector_y,
            unwrapped_phase_vector_z,
            initial_distance_from_oscillation_output_x_fft,
            initial_distance_from_oscillation_output_y_fft,
            initial_distance_from_oscillation_output_z_fft,
            x_fft_initial_y,
            x_fft_initial_z,
            y_fft_initial_y,
            y_fft_initial_z,
            z_fft_initial_y,
            z_fft_initial_z
       )

        push!(simulation_data, data_entry)
    end

    return simulation_data
end