include("helperFunctions.jl")

export
    plotGausWavespeed2d,
    plotGausWavenumber2d,
    plotGausAttenuation2d,
    plotGausAttenuationK2d,
    plot_ωγ_attenuation_2d,
    plot_ωγ_wavespeed_2d,
    plotGuess,
    plotAmp,
    plotPhase,
    plotPhaseRatio,
    plotAmpRatio,
    plotGausPressDep,
    combinePlots,
    plotEllipseAttenuation2d,
    plotAmpTiled,
    combinePlotsTiled,
    getMeanField,
    plotAmpRatioMeanField,
    plotPhaseRatioMeanField,
    plotStitchPhaseScatter,
    plotStitchAttenuation,
    plotStitchAmpRatio,
    plotStitchAmpPhase,
    getMeanField3d


#----------------------------------------------3d specific plots------------------------------------------------------------
function getMeanField3d(filtered_data, transverse_axis; plot = true, shear = false)
    @assert transverse_axis in ["y", "z"] "transverse_axis must be either 'y' or 'z'"
    # fit the prime of the parralell direction data
    if shear == true
        if transverse_axis == "y"
            amplitude_vector = filtered_data[1].amplitude_vector_y
            attenuation = filtered_data[1].alphaoveromega_y
            distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_y_fft
            omega = filtered_data[1].omega # but this is dimensionelss
        else
            amplitude_vector = filtered_data[1].amplitude_vector_z
            attenuation = filtered_data[1].alphaoveromega_z
            distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_z_fft
            omega = filtered_data[1].omega # but this is dimensionelss
        end
    else
        # for compressional waves
        amplitude_vector = filtered_data[1].amplitude_vector_x
        attenuation = filtered_data[1].alphaoveromega_x
        distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_x_fft
        omega = filtered_data[1].omega # but this is dimensionelss
    end

    A = filtered_data[1].pressure/100
    mean_field_amp = A*exp.(-attenuation*omega*distance_from_wall)

    #  This part picks out the "paralell" data
    if shear == true
        if transverse_axis == "y"
            y = filtered_data[1].amplitude_vector_y
            x_parra = filtered_data[1].initial_distance_from_oscillation_output_y_fft
            prime_field_amp = abs.(y - mean_field_amp)
        else
            y = filtered_data[1].amplitude_vector_z
            x_parra = filtered_data[1].initial_distance_from_oscillation_output_z_fft
            prime_field_amp = abs.(y - mean_field_amp)
        end
    else
        y = filtered_data[1].amplitude_vector_x
        x_parra = filtered_data[1].initial_distance_from_oscillation_output_x_fft 
        prime_field_amp = abs.(y - mean_field_amp)
    end

    if plot==true
        mat"""
        coefficients = polyfit($(distance_from_wall), log(abs($(y))), 1);
        fitted_attenuation = coefficients(1);
        intercept_attenuation = coefficients(2);
        mean_field_new = exp(intercept_attenuation) * exp(fitted_attenuation .* $(distance_from_wall));
        prime_field_amp_new = abs($(y)-mean_field_new);
        figure
        scatter($(x_parra), prime_field_amp_new, "*", "DisplayName", " \$ A_{||}' \$")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$A(x)\$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show", "Interpreter", "latex")
        """
    end

    if plot==true
        mat"""
        scatter($(x_parra), $(y), "o", "DisplayName", "\$ A_{||} \$")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$A(x)\$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end
    
    if plot==true
        mat"""
        scatter($(x_parra), mean_field_new, "v","DisplayName", "\$ \\overline{A}_{||} \$")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$A(x)\$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end

    # perpidicular data
    if shear ==true
        x_perp = filtered_data[1].initial_distance_from_oscillation_output_x_fft
        y = filtered_data[1].amplitude_vector_x
    else
        if transverse_axis == "y"
            # !! these are the same length but negative
            x_perp = filtered_data[1].initial_distance_from_oscillation_output_y_fft
            y = filtered_data[1].amplitude_vector_y
            println("x_perp max Amp: ", maximum(x_perp))
            println("y max Amp: ", maximum(y))
        else
            x_perp = filtered_data[1].initial_distance_from_oscillation_output_z_fft
            y = filtered_data[1].amplitude_vector_z
            println("x_perp max Amp: ", maximum(x_perp))
            println("y max Amp: ", maximum(y))
        end
        # x_perp = filtered_data[1].initial_distance_from_oscillation_output_y_fft
        # y = filtered_data[1].amplitude_vector_y
    end

    if plot == true
        mat"""
        scatter($(x_perp), $(y), "o", "DisplayName", "\$ A_\\perp \$")
        set(gca, 'YScale', 'log')
        grid on
        legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
        set(gca, 'FontSize', 15);
        legend('FontSize', 15)
        """
    end

    # phase
    

    # ----------  phase in the x-direction ----------------------
    if shear == true
        phase = filtered_data[1].unwrapped_phase_vector_y 
    else
        phase = filtered_data[1].unwrapped_phase_vector_x
    end
    phase = mod.(phase, 2π)

    if shear == true
        distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_y_fft
    else
        distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    end

    if shear== true
        # wavespeed = filtered_data[1].wavespeed_y # output of sims, but not in crunch
        wavespeed = -omega * .034  /(filtered_data[1].wavenumber_y)# using c = omega_hat * sqrt(K/M) / wavenumber; driving_frequency*2*pi*sqrt(M/K)/(wavenumber*1);
        #  the .03 is a fit factor for the phase
    else
        wavespeed = filtered_data[1].wavespeed_x
    end
    # wavespeed = filtered_data[1].wavespeed_x
    omega = filtered_data[1].omega # but this is dimensionelss
    # wavespeed_x = driving_frequency*2*pi*sqrt(M/K)/(wavenumber*1);

    if shear == true
        wavenumber = filtered_data[1].wavenumber_y
    else
        wavenumber = filtered_data[1].wavenumber_x
    end
    # wavenumber = filtered_data[1].wavenumber_x
    # mean_field = (wavenumber).*distance_from_wall
    
    # _----------------------------- new wrapped feature
    distance_from_wall_sorted, unwrapped_phase_sorted = unwrapScattered(distance_from_wall, phase)
    # create a fit line to distance distance_from_wall_sorted and unwrapped_phase_sorted that are y values that correspond to distance_from_wall_sorted
    fitline = Polynomials.fit(distance_from_wall_sorted, unwrapped_phase_sorted, 1)
    # get the phase of the fit line
    mean_field_phase = fitline.(distance_from_wall_sorted)
    prime_field_phase = unwrapped_phase_sorted .- mean_field_phase
    # -----------------------------

    

    # mean_field_phase = -(omega/wavespeed).*distance_from_wall
    mean_field_phase = mod.(mean_field_phase, 2π)
    # prime_field_phase = phase .- mean_field_phase
    prime_field_phase = mod.(prime_field_phase, 2π)
    prime_field_amp_new = prime_field_phase.- mean(prime_field_phase);
    
    if shear == true
        prime_field_amp_new = abs.(prime_field_amp_new) # Did this because it was negative
        # prime_field_amp_new = prime_field_phase
    end
    if plot==true
        mat"""
        figure
        scatter($(distance_from_wall), $(prime_field_amp), "*", "DisplayName", "\$ \\phi_{||}' \$")
        hold on
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$ \\phi(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end


    if plot==true
        mat"""
        scatter($(distance_from_wall), $(phase), "o", "DisplayName", "\$ \\phi_{||} \$")
        hold on
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$\\phi(x \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show", "interpreter", "latex")
        """
    end
    
    if plot==true
        mat"""
        scatter($(distance_from_wall_sorted), $(mean_field_phase), "v","DisplayName", "\$ \\overline{\\phi}_{||} \$ ")
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$ \\phi(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show", "interpreter", "latex")
        """
    end
   


# ---------------------------------- this is all for the transverse data--------------------------------
   if shear == true
        x_perp = filtered_data[1].initial_distance_from_oscillation_output_x_fft
        y = filtered_data[1].unwrapped_phase_vector_x
    else
        if transverse_axis == "y"
            # !!! Problem, these are differen lengths !!!
            x_perp = filtered_data[1].initial_distance_from_oscillation_output_y_fft
            y = filtered_data[1].unwrapped_phase_vector_y 
            println("length of x_perp: ", length(x_perp))
            println("length of y:", length(y))
        else
            x_perp = filtered_data[1].initial_distance_from_oscillation_output_z_fft
            y = filtered_data[1].unwrapped_phase_vector_z
            println("length of x_perp: ", length(x_perp))
            println("length of y:", length(y))
        end
    end 

    y = mod.(y, 2π)
    
    
    if plot == true
        mat"""
        scatter($(x_perp), $(y), "o", "DisplayName", " \$ \\phi_\\perp \$ ")
        grid on
        legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
        set(gca, 'FontSize', 15)
        legend('FontSize', 15, "interpreter", "latex")
        """
    end

    # new_y =  y + filtered_data[1].wavenumber_x.*x_perp
    # new_y = mod.(new_y, 2π)
    # if plot == true
    #     mat"""
    #     scatter($(x_perp), $(new_y), "o", "DisplayName", "raw -  mean-phase")
    #     grid on
    #     legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
    #     set(gca, 'FontSize', 15)
    #     legend('FontSize', 15)
    #     """
    # end

    return mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase
end
#-----------------------------------------------Stitched PLots------------------------------------------------------------
function plotStitchAmpPhase(simulation_data, gamma_values) 
    mat"""
    ax_energy = figure;
    %xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\hat{\\gamma} \\frac{2\\hat{c}^2 \\left( 1 - \\cos\\overline{\\sigma}_{\\Delta \\phi_{ij}} \\right)}{\\hat{d^2\\omega}^2}\\left(\\overline{\\frac{A_{\\perp}}{A_{\\parallel}}}\\right)^2 \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;
    """
    marker_shape_vector = ["-*", "-o", "-v", "-+", "-.", "-x", "d"]

    for γ_value in gamma_values
        # filter the data based on those that are close to gamma_value
        closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
        closest_γ_value = simulation_data[closest_γ_index].gamma
        matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
        plot_gamma = γ_value
        gamma_value = γ_value
        max_gamma = maximum(gamma_values)

        # Get a list of unique input pressures
        pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry
        # closest_p_idx = argmin(abs.(pressure_list .- .05))
        # pressure_list = [minimum(pressure_list), pressure_list[closest_p_idx], maximum(pressure_list)] # just get the limits
        pressure_list = [minimum(pressure_list), maximum(pressure_list)] # just get the limits

        # get a range for plotting color from 0 to 1
        normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

        # Create a line for each pressure
        for pressure_value in pressure_list

            # Assign a color
            idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
            marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

            # Only look at data for current pressure value
            matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Initizalized vectors for just this pressure
            loop_mean_E_list = Float64[];
            loop_mean_attenuation_list = Float64[];

            # Look at a single omega gamma value since each one spans all seeds
            matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
            
            for omega_gamma_value in matching_omega_gamma_list

                # Only look at data for current omega_gamma value
                matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

                # Get the mean over all seeds       
                jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
                E_ratio_list = Float64[]
                seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))

                for k_seed in seed_list
                    k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                    k_seed_omega = k_seed_data[1].omega
                    phase_vector_y = k_seed_data[1].unwrapped_phase_vector_y
                    # Wrap the phase vector around 2π
                    wrapped_phase = mod.(phase_vector_y, 2π)
                    distance_from_wall = k_seed_data[1].initial_distance_from_oscillation_output_y_fft
                    # mean_distance = meanDistNeighbor(distance_from_wall, wrapped_phase)
                    if isempty(k_seed_data[1].unwrapped_phase_vector_y)
                        println("Empty y-phase vector for: Pressure $(pressure_value) OmegaGamma $(omega_gamma_value) seed $(k_seed)")
                        continue
                    end
                    mean_distance = plotAmp(k_seed_data; plot=false)
                    mean_phase = plotPhase(k_seed_data; plot=false)
                    mean_scatter = 1-cos(mean_phase)
                    # mean_distance = mean_scatter.*mean_distance.^2
                    mean_distance = mean_scatter.*mean_distance.^2 ./ k_seed_omega^2 ./ 1.2^2 .* 2
                    push!(E_ratio_list, mean_distance)
                end

                j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
                push!(loop_mean_E_list, j_E_ratio)
                push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
            end

            # This is needed because MATLAB.jl has a hard time escaping \'s
            pressure_label = @sprintf("\$ %.4f, %.4f \$", gamma_value, pressure_value)

            gamma_val = γ_value
            marker_shape = marker_shape_vector[findfirst(==(gamma_val), gamma_values)]
            mat"""
            omega_gamma = $(matching_omega_gamma_list);
            loop_mean_E_list = $(loop_mean_E_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            plot_gamma = $(plot_gamma);
            marker_color = $(marker_color);
            pressure_label = $(pressure_label);
            marker_shape = $(marker_shape)
            marker_size = exp(plot_gamma/$(max_gamma))*3

            % plot( omega_gamma/$(gamma_val), loop_mean_E_list, marker_shape, 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            % plot( omega_gamma, loop_mean_E_list, marker_shape, 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            plot( omega_gamma, loop_mean_E_list, '-o', 'MarkerSize', marker_size,'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            % plot( omega_gamma/plot_gamma, loop_mean_E_list, '-o', 'MarkerSize', marker_size,'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            """
        end

    end
    # Add legends to the plots
    mat"""
    % legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    title(leg, "\$  \\hat{\\gamma}, \\hat{P} \$")
    fitx = [.003, 2]
    fity = .001*fitx.^-(2/3)
    fitz = .03*fitx.^-(2/3)
    leg.AutoUpdate = 'off'; 
    plot(fitx, fity, 'k-', 'LineWidth', 3, 'DisplayName', '');  % No legend for fity
    plot(fitx, fitz, 'k-', 'LineWidth', 3, 'DisplayName', '');  % No legend for fitz
    text(.003, .4, '\$ -\\frac{2}{3} \$', 'Interpreter', 'latex', 'FontSize', 20);
    text(.003, .01, '\$ -\\frac{2}{3} \$', 'Interpreter', 'latex', 'FontSize', 20);
    """ 
end
function plotStitchAmpRatio(simulation_data, gamma_values; shear=false) 
    mat"""
    ax_energy = figure;
    xlabel('\$  \\hat{\\omega} \$', "FontSize", 20, "Interpreter", "latex");

    %xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\hat{\\gamma} \\left( \\overline{\\frac{A_{\\perp}}{A_{\\parallel}}} \\right) ^2\$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;
    """
    marker_shape_vector = ["-*", "-o", "-v", "-+", "-.", "-x", "d"]

    for γ_value in gamma_values
        # filter the data based on those that are close to gamma_value
        closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
        closest_γ_value = simulation_data[closest_γ_index].gamma
        matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
        plot_gamma = γ_value
        gamma_value = γ_value
        max_gamma = maximum(gamma_values)

        # Get a list of unique input pressures
        pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry
        pressure_list = [minimum(pressure_list), maximum(pressure_list)] # just get the limits

        # get a range for plotting color from 0 to 1
        normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

        # Create a line for each pressure
        for pressure_value in pressure_list

            # Assign a color
            idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
            marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

            # Only look at data for current pressure value
            matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Initizalized vectors for just this pressure
            loop_mean_E_list = Float64[];
            loop_mean_attenuation_list = Float64[];

            # Look at a single omega gamma value since each one spans all seeds
            matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
            
            for omega_gamma_value in matching_omega_gamma_list

                # Only look at data for current omega_gamma value
                matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

                # Get the mean over all seeds       
                jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
                E_ratio_list = Float64[]
                seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))

                for k_seed in seed_list
                    k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                    k_seed_omega = k_seed_data[1].omega
                    phase_vector_y = k_seed_data[1].unwrapped_phase_vector_y
                    # Wrap the phase vector around 2π
                    wrapped_phase = mod.(phase_vector_y, 2π)
                    distance_from_wall = k_seed_data[1].initial_distance_from_oscillation_output_y_fft
                    # mean_distance = meanDistNeighbor(distance_from_wall, wrapped_phase)
                    if isempty(k_seed_data[1].unwrapped_phase_vector_y)
                        println("Empty y-phase vector for: Pressure $(pressure_value) OmegaGamma $(omega_gamma_value) seed $(k_seed)")
                        continue
                    end
                    mean_distance = plotAmp(k_seed_data; plot=false, shear=shear)
                    push!(E_ratio_list, mean_distance)
                end

                j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
                push!(loop_mean_E_list, j_E_ratio)
                push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
            end

            # This is needed because MATLAB.jl has a hard time escaping \'s
            pressure_label = @sprintf("\$ %.4f, %.4f \$", pressure_value, gamma_value)

            gamma_val = γ_value
            marker_shape = marker_shape_vector[findfirst(==(gamma_val), gamma_values)]
            mat"""
            omega_gamma = $(matching_omega_gamma_list);
            loop_mean_E_list = $(loop_mean_E_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            plot_gamma = $(plot_gamma);
            marker_color = $(marker_color);
            pressure_label = $(pressure_label);
            marker_shape = $(marker_shape);
            marker_size = exp(plot_gamma/$(max_gamma))*3
            y =  plot_gamma*loop_mean_E_list.^2;
            y(y > 1) = NaN;
            plot( omega_gamma/$(gamma_val),y, "-o", 'MarkerSize', marker_size , 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            %plot( omega_gamma, loop_mean_E_list, marker_shape, 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            """
        end

    end
    # Add legends to the plots
    mat"""
    % legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    title(leg, "\$  \\hat{P}, \\hat{\\gamma} \$")
    fitx = [.03, 2]
    fity = .2*fitx.^1
    fitz = .001*fitx.^1
    leg.AutoUpdate = 'off'; 
    plot(fitx, fity, 'k-', 'LineWidth', 3, 'DisplayName', 'slope = 1')
    plot(fitx, fitz, 'k-', 'LineWidth', 3, 'DisplayName', 'slope = 1')
    set(get(gca, 'ylabel'), 'rotation', 0);
    text(.2, .000, '\$ 1 \$', 'Interpreter', 'latex', 'FontSize', 20);
    text(.2, .08, '\$ 1 \$', 'Interpreter', 'latex', 'FontSize', 20);
    """ 
end

function plotStitchAttenuation(simulation_data, gamma_values, mean_diameter; shear=false) 
    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    mat"""
    figure_attenuation = figure;
    loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    """
    marker_shape_vector = ["-*", "-o", "-v", "-+", "-.", "-x", "-d"]
    if shave
        marker_shape_vector = ["*", "o", "v", "+", ".", "x", "d"]
    end
    for γ_value in gamma_values
        # filter the data based on those that are close to gamma_value
        closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
        closest_γ_value = simulation_data[closest_γ_index].gamma
        matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
        plot_gamma = γ_value
        gamma_value = γ_value

        # Get a list of unique input pressures
        pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry
        pressure_list = [minimum(pressure_list), maximum(pressure_list)] # just get the limits

        # get a range for plotting color from 0 to 1
        normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

        # Create a line for each pressure
        for pressure_value in pressure_list

            # Assign a color
            idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
            marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

            # Only look at data for current pressure value
            matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Initizalized vectors for just this pressure
            loop_mean_E_list = Float64[];
            loop_mean_attenuation_list = Float64[];

            # Look at a single omega gamma value since each one spans all seeds
            matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
            
            for omega_gamma_value in matching_omega_gamma_list

                # Only look at data for current pressure value
                matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
    
                # Get the mean over all seeds
                loop_mean_alphaoveromega = mean_diameter .* mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
    
                # Append values
                push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
            end

            # This is needed because MATLAB.jl has a hard time escaping \'s
            pressure_label = @sprintf("\$ %.4f, %.4f \$", pressure_value, gamma_value)

            gamma_val = γ_value
            marker_shape = marker_shape_vector[findfirst(==(gamma_val), gamma_values)]
            if shave
                middle_shave_length = 2
                start_index = div(length(loop_mean_attenuation_list) - middle_shave_length, 2) + 1
                loop_mean_attenuation_list = loop_mean_attenuation_list[start_index:start_index + (middle_shave_length - 1)]
                matching_omega_gamma_list = matching_omega_gamma_list[start_index:start_index + (middle_shave_length - 1)]
            end
            mat"""
            x = $(matching_omega_gamma_list);
            y = $(loop_mean_attenuation_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            plot_gamma = $(plot_gamma);
            marker_color = $(marker_color);
            pressure_label = $(pressure_label);
            marker_shape = $(marker_shape);
            
            plot( x, y, marker_shape, 'Color', marker_color, 'DisplayName', pressure_label);
            %plot( omega_gamma, loop_mean_E_list, marker_shape, 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            """
        end


    end
    # Add legends to the plots
    mat"""
    % legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    title(leg, "\$  \\hat{P}, \\hat{\\gamma} \$")
    """ 
end
function plotStitchPhaseScatter(simulation_data, gamma_values; shear=false) 
    mat"""
    ax_energy = figure;
    xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
    %xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ 1-\\cos \\overline{\\sigma}_{\\Delta \\phi_{\\perp}} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;
    """
    marker_shape_vector = ["-*", "-o", "-v", "-+", "-.", "-x", "d"]

    for γ_value in gamma_values
        # filter the data based on those that are close to gamma_value
        closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
        closest_γ_value = simulation_data[closest_γ_index].gamma
        matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
        plot_gamma = γ_value
        gamma_value = γ_value
        max_gamma = maximum(gamma_values)

        # Get a list of unique input pressures
        pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry
        pressure_list = [minimum(pressure_list), maximum(pressure_list)] # just get the limits

        # get a range for plotting color from 0 to 1
        normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

        # Create a line for each pressure
        for pressure_value in pressure_list

            # Assign a color
            idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
            marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

            # Only look at data for current pressure value
            matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Initizalized vectors for just this pressure
            loop_mean_E_list = Float64[];
            loop_mean_attenuation_list = Float64[];

            # Look at a single omega gamma value since each one spans all seeds
            matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
            
            for omega_gamma_value in matching_omega_gamma_list

                # Only look at data for current omega_gamma value
                matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

                # Get the mean over all seeds       
                jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
                E_ratio_list = Float64[]
                seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))

                for k_seed in seed_list
                    k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                    # k_seed_omega = k_seed_data[1].omega
                    # phase_vector_y = k_seed_data[1].unwrapped_phase_vector_y
                    # Wrap the phase vector around 2π
                    # wrapped_phase = mod.(phase_vector_y, 2π)
                    # distance_from_wall = k_seed_data[1].initial_distance_from_oscillation_output_y_fft
                    # mean_distance = meanDistNeighbor(distance_from_wall, wrapped_phase)
                    if isempty(k_seed_data[1].unwrapped_phase_vector_y)
                        println("Empty y-phase vector for: Pressure $(pressure_value) OmegaGamma $(omega_gamma_value) seed $(k_seed)")
                        continue
                    end
                    mean_scatter = plotPhase(k_seed_data; plot=false, shear=shear)
                    mean_scatter = isinf(mean_scatter) ? NaN : mean_scatter # saftey for infinite values
                    push!(E_ratio_list, 1-cos(mean_scatter))
                    # push!(E_ratio_list, mean_distance)
                end

                j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
                push!(loop_mean_E_list, j_E_ratio)
                push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
            end

            # This is needed because MATLAB.jl has a hard time escaping \'s
            pressure_label = @sprintf("\$ %.4f, %.4f \$", pressure_value, gamma_value)

            gamma_val = γ_value
            marker_shape = marker_shape_vector[findfirst(==(gamma_val), gamma_values)]
            mat"""
            omega_gamma = $(matching_omega_gamma_list);
            loop_mean_E_list = $(loop_mean_E_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            plot_gamma = $(plot_gamma);
            marker_color = $(marker_color);
            pressure_label = $(pressure_label);
            marker_shape = $(marker_shape);
            marker_size = exp(plot_gamma/$(max_gamma))*3;

            plot( omega_gamma/$(gamma_val), loop_mean_E_list, "-o", 'MarkerSize', marker_size, 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            %plot( omega_gamma, loop_mean_E_list, marker_shape, 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
            """
        end


    end
    # Add legends to the plots
    mat"""
    % legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    title(leg, "\$  \\hat{P}, \\hat{\\gamma} \$")
    fitx = [.03, .8]
    fity = 2*fitx.^.5
    fitz = 1*fitx.^2
    leg.AutoUpdate = 'off'; 
    plot(fitx, fity, 'k-', 'LineWidth', 3, 'DisplayName', 'slope = 1/2')
    plot(fitx, fitz, 'k-',  'LineWidth', 3,'DisplayName', 'slope = 2')
    text(.1, 1.2, '\$ \\frac{1}{2} \$', 'Interpreter', 'latex', 'FontSize', 20);
    text(.2, .02, '\$ 2 \$', 'Interpreter', 'latex', 'FontSize', 20)
    """ 
end

#-----------------------------------------------Mean Field------------------------------------------------------------
function plotPhaseRatioMeanField(simulation_data, γ_value) 

    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
    plot_gamma = γ_value
    gamma_value = γ_value

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry

    # Limit range to data
    upper_limit_line_x = [1*γ_value; 1*γ_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*γ_value; .1*γ_value]
    lower_limit_line_y = [1E-5; 1]

   mat"""
    ax_energy = figure;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ 1-\\cos \\overline{\\sigma}_{\\Delta \\phi_{\\perp}} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;
    """

    # get a range for plotting color from 0 to 1
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_E_list = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
        
        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current omega_gamma value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds       
            jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
            E_ratio_list = Float64[]
            seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))
            for k_seed in seed_list
                k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                k_seed_omega = k_seed_data[1].omega
                phase_vector_y = k_seed_data[1].unwrapped_phase_vector_y

                # Wrap the phase vector around 2π
                wrapped_phase = mod.(phase_vector_y, 2π)
                distance_from_wall = k_seed_data[1].initial_distance_from_oscillation_output_y_fft
                # mean_distance = meanDistNeighbor(distance_from_wall, wrapped_phase)
                if isempty(k_seed_data[1].unwrapped_phase_vector_y)
                    println("Empty y-phase vector for: Pressure $(pressure_value) OmegaGamma $(omega_gamma_value) seed $(k_seed)")
                    continue
                end
                mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(k_seed_data; plot=false)
                distance_y = k_seed_data[1].initial_distance_from_oscillation_output_x_fft
                scatter_y = meanPhaseDev(distance_y, prime_field_phase, 1)
                mean_distance = scatter_y
                # mean_distance = scatter_y - mean(mean_distance)
                mean_distance = isinf(mean_distance) ? NaN : mean_distance # saftey for infinite values
                push!(E_ratio_list, 1-cos(mean_distance))
            end
            j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
            push!(loop_mean_E_list, j_E_ratio)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end

       @bp
        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$ %.4f\$", pressure_value)

        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        loop_mean_E_list = $(loop_mean_E_list);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);

        % Plot Attenuation
        %loglog(ax_attenuation, omega_gamma, mean_attenuation_x, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
        % Plot Aspect Ratio
        % plot(ax_energy, omega_gamma, loop_mean_E_list, 'o-', 'Color', marker_color, 'DisplayName', pressure_label);
        % plot( omega_gamma, loop_mean_E_list, 'o-', 'Color', marker_color, 'DisplayName', pressure_label);
        plot( omega_gamma, loop_mean_E_list, 'o-', 'Color', marker_color);
        """
    end

    # Add legends to the plots
    mat"""
    % legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    % leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    % title(leg, "\$ \\hat{P} \$")
    """
end
function plotAmpRatioMeanField(simulation_data, γ_value)

    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
    plot_gamma = γ_value
    gamma_value = γ_value

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry

    # Limit range to data
    upper_limit_line_x = [1*γ_value; 1*γ_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*γ_value; .1*γ_value]
    lower_limit_line_y = [1E-5; 1]

    # # Start MATLAB session
    # mat"""
    # figure_main = figure;
    # %tiled_main = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'none'); % 2 rows, 1 column

    # % Axes for Attenuation
    # ax_attenuation = nexttile;
    # hold(ax_attenuation, 'on');
    # % xlabel(ax_attenuation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    # ylabel(ax_attenuation, '\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}}\$', "FontSize", 20, "Interpreter", "latex");
    # set(ax_attenuation, 'XScale', 'log');
    # set(ax_attenuation, 'YScale', 'log')
    # set(get(ax_attenuation, 'ylabel'), 'rotation', 0);
    # grid(ax_attenuation, 'on');
    # box(ax_attenuation, 'on');
    # %plot(ax_attenuation, $(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$');
    # %plot(ax_attenuation, $(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$');
    # set(ax_attenuation, 'XTickLabel', []);

    # % Axes for Energy
    # ax_energy = nexttile;
    # hold(ax_energy, 'on');
    # xlabel(ax_energy, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    # ylabel(ax_energy, '\$ \\overline{\\frac{A_{\\parallel}}{A_{\\perp}}} \$', "FontSize", 15, "Interpreter", "latex");
    # set(ax_energy, 'XScale', 'log');
    # set(ax_energy, 'YScale', 'log')
    # set(get(ax_energy, 'ylabel'), 'rotation', 0);
    # grid(ax_energy, 'on');
    # box(ax_energy, 'on');
    # %set(ax_energy, 'XTickLabel', []);
    # """

    mat"""
    ax_energy = figure;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$  \\overline{\\frac{A_{\\perp}}{A_{\\parallel}}} \$', "FontSize", 20, "Interpreter", "latex");
    set(get(gca, 'ylabel'), 'rotation', 0);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;
    """

    # get a range for plotting color from 0 to 1
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_E_list = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
        
        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current omega_gamma value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds       
            jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
            E_ratio_list = Float64[]
            seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))
            for k_seed in seed_list
                k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                k_seed_omega = k_seed_data[1].omega
                phase_vector_y = k_seed_data[1].unwrapped_phase_vector_y
                # Wrap the phase vector around 2π
                wrapped_phase = mod.(phase_vector_y, 2π)
                distance_from_wall = k_seed_data[1].initial_distance_from_oscillation_output_y_fft
                # mean_distance = meanDistNeighbor(distance_from_wall, wrapped_phase)
                mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase = getMeanField(k_seed_data; plot=false)
                mean_distance = mean(prime_field_amp ./ mean_field_amp)
                push!(E_ratio_list, mean_distance)
            end
            j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
            push!(loop_mean_E_list, j_E_ratio)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        # valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        # matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        # loop_mean_E_list = loop_mean_E_list[valid_indices]
        # loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]
        @bp
        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$ %.4f \$", pressure_value)

        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        loop_mean_E_list = $(loop_mean_E_list);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);

        % loglog( omega_gamma, loop_mean_E_list, 'o-', 'Color', marker_color, 'DisplayName', pressure_label);
        loglog( omega_gamma, loop_mean_E_list, 'o-', 'Color', marker_color);
        ylim([.01, .6])
        xlim([min(omega_gamma), 1])
        """
    end

    # Add legends to the plots
    mat"""
    % leg = legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
    % title(leg, "\$ \\hat{P} \$")
    """ 
end
    
function plotAmpRatio(simulation_data, γ_value)
   
    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
    plot_gamma = γ_value
    gamma_value = γ_value

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry

    # Limit range to data
    upper_limit_line_x = [1*γ_value; 1*γ_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*γ_value; .1*γ_value]
    lower_limit_line_y = [1E-5; 1]

    # # Start MATLAB session
    # mat"""
    # figure_main = figure;
    # %tiled_main = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'none'); % 2 rows, 1 column

    # % Axes for Attenuation
    # ax_attenuation = nexttile;
    # hold(ax_attenuation, 'on');
    # % xlabel(ax_attenuation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    # ylabel(ax_attenuation, '\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}}\$', "FontSize", 20, "Interpreter", "latex");
    # set(ax_attenuation, 'XScale', 'log');
    # set(ax_attenuation, 'YScale', 'log')
    # set(get(ax_attenuation, 'ylabel'), 'rotation', 0);
    # grid(ax_attenuation, 'on');
    # box(ax_attenuation, 'on');
    # %plot(ax_attenuation, $(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$');
    # %plot(ax_attenuation, $(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$');
    # set(ax_attenuation, 'XTickLabel', []);

    # % Axes for Energy
    # ax_energy = nexttile;
    # hold(ax_energy, 'on');
    # xlabel(ax_energy, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    # ylabel(ax_energy, '\$ \\overline{\\frac{A_{\\parallel}}{A_{\\perp}}} \$', "FontSize", 15, "Interpreter", "latex");
    # set(ax_energy, 'XScale', 'log');
    # set(ax_energy, 'YScale', 'log')
    # set(get(ax_energy, 'ylabel'), 'rotation', 0);
    # grid(ax_energy, 'on');
    # box(ax_energy, 'on');
    # %set(ax_energy, 'XTickLabel', []);
    # """

    mat"""
    ax_energy = figure;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$  \\overline{\\frac{A_{\\perp}}{A_{\\parallel}}} \$', "FontSize", 20, "Interpreter", "latex");
    set(get(gca, 'ylabel'), 'rotation', 0);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;
    """

    # get a range for plotting color from 0 to 1
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_E_list = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
        
        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current omega_gamma value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds       
            jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
            E_ratio_list = Float64[]
            seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))
            for k_seed in seed_list
                k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                k_seed_omega = k_seed_data[1].omega
                phase_vector_y = k_seed_data[1].unwrapped_phase_vector_y
                # Wrap the phase vector around 2π
                wrapped_phase = mod.(phase_vector_y, 2π)
                distance_from_wall = k_seed_data[1].initial_distance_from_oscillation_output_y_fft
                # mean_distance = meanDistNeighbor(distance_from_wall, wrapped_phase)
                mean_distance = plotAmp(k_seed_data; plot=false)
                push!(E_ratio_list, mean_distance)
            end
            j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
            push!(loop_mean_E_list, j_E_ratio)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        # valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        # matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        # loop_mean_E_list = loop_mean_E_list[valid_indices]
        # loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]
        @bp
        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$\\hat{P} = %.4f\$", pressure_value)

        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        loop_mean_E_list = $(loop_mean_E_list);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);

        % Plot Attenuation
        %loglog(ax_attenuation, omega_gamma, mean_attenuation_x, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
        % Plot Aspect Ratio
        loglog( omega_gamma, loop_mean_E_list, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        ylim([.01, .6])
        xlim([min(omega_gamma), 1])
        """
    end

    # Add legends to the plots
    mat"""
    %legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
    """ 
        
end


function getMeanField(filtered_data; plot = true, shear = false)
    if shear == true
        amplitude_vector = filtered_data[1].amplitude_vector_y
        attenuation = filtered_data[1].alphaoveromega_y
        distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_y_fft
        omega = filtered_data[1].omega # but this is dimensionelss
    else
        amplitude_vector = filtered_data[1].amplitude_vector_x
        attenuation = filtered_data[1].alphaoveromega_x
        distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_x_fft
        omega = filtered_data[1].omega # but this is dimensionelss
    end
    # this was the old version before adding shear
    # attenuation = filtered_data[1].alphaoveromega_x # will need to figure out how to adapt this for y later
    # distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    # omega = filtered_data[1].omega # but this is dimensionelss

    A = filtered_data[1].pressure/100
    mean_field_amp = A*exp.(-attenuation*omega*distance_from_wall)
    # log_amplitude = log.(abs(amplitude_vector_x))
    # coefficents = polyfit(distance_from_wall, log_amplitude, 1)
    if shear == true
        y = filtered_data[1].amplitude_vector_y
        x_parra = filtered_data[1].initial_distance_from_oscillation_output_y_fft 
        prime_field_amp = abs.(y - mean_field_amp)
    else
        y = filtered_data[1].amplitude_vector_x
        x_parra = filtered_data[1].initial_distance_from_oscillation_output_x_fft 
        prime_field_amp = abs.(y - mean_field_amp)
    end
    # old version of above 
    # y = filtered_data[1].amplitude_vector_x    
    # prime_field_amp = abs.(y - mean_field_amp)

    if plot==true
        mat"""
        coefficients = polyfit($(distance_from_wall), log(abs($(y))), 1);
        fitted_attenuation = coefficients(1);
        intercept_attenuation = coefficients(2);
        mean_field_new = exp(intercept_attenuation) * exp(fitted_attenuation .* $(distance_from_wall));
        prime_field_amp_new = abs($(y)-mean_field_new);
        figure
        scatter($(x_parra), prime_field_amp_new, "*", "DisplayName", " \$ A_{||}' \$")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$A(x)\$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show", "Interpreter", "latex")
        """
    end


    if plot==true
        mat"""
        scatter($(x_parra), $(y), "o", "DisplayName", "\$ A_{||} \$")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$A(x)\$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end
    
    if plot==true
        mat"""
        scatter($(x_parra), mean_field_new, "v","DisplayName", "\$ \\overline{A}_{||} \$")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$A(x)\$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end
    
    if shear ==true
        x_perp = filtered_data[1].initial_distance_from_oscillation_output_x_fft
        y = filtered_data[1].amplitude_vector_x
    else
        x_perp = filtered_data[1].initial_distance_from_oscillation_output_y_fft
        y = filtered_data[1].amplitude_vector_y
    end

    if plot == true
        mat"""
        scatter($(x_perp), $(y), "o", "DisplayName", "\$ A_\\perp \$")
        set(gca, 'YScale', 'log')
        grid on
        legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
        set(gca, 'FontSize', 15);
        legend('FontSize', 15)
        """
    end

    # phase
    if shear == true
        phase = filtered_data[1].unwrapped_phase_vector_y 
    else
        phase = filtered_data[1].unwrapped_phase_vector_x
    end
    # phase = filtered_data[1].unwrapped_phase_vector_x # will need to figure out how to adapt this for y later
    phase = mod.(phase, 2π)
    if shear == true
        distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_y_fft
    else
        distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    end
    # distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    if shear== true
        # wavespeed = filtered_data[1].wavespeed_y # output of sims, but not in crunch
        wavespeed = -omega * .034  /(filtered_data[1].wavenumber_y)# using c = omega_hat * sqrt(K/M) / wavenumber; driving_frequency*2*pi*sqrt(M/K)/(wavenumber*1);
        #  the .03 is a fit factor for the phase
    else
        wavespeed = filtered_data[1].wavespeed_x
    end
    # wavespeed = filtered_data[1].wavespeed_x
    omega = filtered_data[1].omega # but this is dimensionelss
    # wavespeed_x = driving_frequency*2*pi*sqrt(M/K)/(wavenumber*1);
    if shear == true
        wavenumber = filtered_data[1].wavenumber_y
    else
        wavenumber = filtered_data[1].wavenumber_x
    end
    # wavenumber = filtered_data[1].wavenumber_x
    # mean_field = (wavenumber).*distance_from_wall
    
    # _-----------------------------
    distance_from_wall_sorted, unwrapped_phase_sorted = unwrapScattered(distance_from_wall, phase)
    # create a fit line to distance distance_from_wall_sorted and unwrapped_phase_sorted that are y values that correspond to distance_from_wall_sorted
    fitline = Polynomials.fit(distance_from_wall_sorted, unwrapped_phase_sorted, 1)
    # get the phase of the fit line
    mean_field_phase = fitline.(distance_from_wall_sorted)
    prime_field_phase = unwrapped_phase_sorted .- mean_field_phase
    # _-----------------------------



    # mean_field_phase = -(omega/wavespeed).*distance_from_wall
    mean_field_phase = mod.(mean_field_phase, 2π)
    # prime_field_phase = phase .- mean_field_phase
    prime_field_phase = mod.(prime_field_phase, 2π)
    prime_field_amp_new = prime_field_phase.- mean(prime_field_phase);
    
    if shear == true
        prime_field_amp_new = abs.(prime_field_amp_new) # Did this because it was negative
        # prime_field_amp_new = prime_field_phase
    end
    if plot==true
        mat"""
        figure
        scatter($(distance_from_wall), $(prime_field_amp), "*", "DisplayName", "\$ \\phi_{||}' \$")
        hold on
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$ \\phi(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end


    if plot==true
        mat"""
        scatter($(distance_from_wall), $(phase), "o", "DisplayName", "\$ \\phi_{||} \$")
        hold on
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$\\phi(x \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show", "interpreter", "latex")
        """
    end
    
    if plot==true
        mat"""
        scatter($(distance_from_wall_sorted), $(mean_field_phase), "v","DisplayName", "\$ \\overline{\\phi}_{||} \$ ")
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$ \\phi(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show", "interpreter", "latex")
        """
    end
   if shear == true
        x_perp = filtered_data[1].initial_distance_from_oscillation_output_x_fft
        y = filtered_data[1].unwrapped_phase_vector_x
    else
        x_perp = filtered_data[1].initial_distance_from_oscillation_output_y_fft
        y = filtered_data[1].unwrapped_phase_vector_y
    end 
    # x_perp = filtered_data[1].initial_distance_from_oscillation_output_y_fft
    if shear == true
        y = filtered_data[1].unwrapped_phase_vector_x
    else
        y = filtered_data[1].unwrapped_phase_vector_y
    end
    # y = filtered_data[1].unwrapped_phase_vector_y
    y = mod.(y, 2π)
    
    transverse_fitline = -1/6 .* x_perp # this is for the lower pressure data = FilterData(simulation_data, .001, :pressure, .1, :omega, .5, :gamma, 1, :seed)
    # transverse_fitline = 1/7 .* x_perp # for the higher pressure     data = FilterData(simulation_data, .1, :pressure, .1, :omega, .5, :gamma, 1, :seed)

    transverse_fitline = mod.(transverse_fitline, 2π)
    z = abs.(transverse_fitline .- y)

    
    if plot == true
        mat"""
        scatter($(x_perp), $(y), "o", "DisplayName", " \$ \\phi_\\perp \$ ")
        grid on
        legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
        set(gca, 'FontSize', 15)
        % plot($(x_perp),$(transverse_fitline), "o", "DisplayName", "Transverse Fit")
        % plot($(x_perp), $(z), "o", "DisplayName", "Transverse Fit-Data")
        legend('FontSize', 15, "interpreter", "latex")
        """
    end

    # new_y =  y + filtered_data[1].wavenumber_x.*x_perp
    # new_y = mod.(new_y, 2π)
    # if plot == true
    #     mat"""
    #     scatter($(x_perp), $(new_y), "o", "DisplayName", "raw -  mean-phase")
    #     grid on
    #     legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
    #     set(gca, 'FontSize', 15)
    #     legend('FontSize', 15)
    #     """
    # end

    return mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase
end

function combinePlotsTiled()
    mat"combinePlotsTiled('fig1.fig', 'fig2.fig');"
end
#-----------------------------------Ellipse Plots------------------------------------------------------------------------
function plotEllipseAttenuation2d(simulation_data, γ_value) 

    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
    plot_gamma = γ_value
    gamma_value = γ_value

    # Get a list of unique input pressures
    pressure_list = unique([entry.pressure for entry in matching_γ_data]) # goes through each entry of simulation_data and get the P value at that entry

    # Limit range to data
    upper_limit_line_x = [1*γ_value; 1*γ_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*γ_value; .1*γ_value]
    lower_limit_line_y = [1E-5; 1]

    # Start MATLAB session
    mat"""
    figure_main = figure;
    tiled_main = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'none'); % 3 rows, 1 column

    % Axes for Attenuation
    ax_attenuation = nexttile;
    hold(ax_attenuation, 'on');
    % xlabel(ax_attenuation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel(ax_attenuation, '\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}}\$', "FontSize", 20, "Interpreter", "latex");
    set(ax_attenuation, 'XScale', 'log');
    set(ax_attenuation, 'YScale', 'log')
    set(get(ax_attenuation, 'ylabel'), 'rotation', 0);
    grid(ax_attenuation, 'on');
    box(ax_attenuation, 'on');
    %plot(ax_attenuation, $(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$');
    %plot(ax_attenuation, $(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$');
    set(ax_attenuation, 'XTickLabel', []);

    % Axes for Rotation Angle
    ax_rotation = nexttile;
    hold(ax_rotation, 'on');
    % xlabel(ax_rotation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel(ax_rotation, '\$ \\overline{ \\left| \\theta \\right| }  \$', "FontSize", 20, "Interpreter", "latex");
    set(ax_rotation, 'XScale', 'log');
    set(get(ax_rotation, 'ylabel'), 'rotation', 0);
    grid(ax_rotation, 'on');
    box(ax_rotation, 'on');
    set(ax_rotation, 'XTickLabel', []);

    % Axes for Aspect Ratio
    ax_aspect_ratio = nexttile;
    hold(ax_aspect_ratio, 'on');
    xlabel(ax_aspect_ratio, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel(ax_aspect_ratio, '\$ \\overline{\\frac{b}{a}} \$', "FontSize", 20, "Interpreter", "latex");
    set(ax_aspect_ratio, 'XScale', 'log');
    set(get(ax_aspect_ratio, 'ylabel'), 'rotation', 0);
    grid(ax_aspect_ratio, 'on');
    box(ax_aspect_ratio, 'on');
    """

    # get a range for plotting color from 0 to 1
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_aspect_ratio_list = Float64[];
        loop_mean_rotation_angles = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))

        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current pressure value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            jvalue_mean_aspect_ratio = mean(entry.mean_aspect_ratio for entry in matching_omega_gamma_data) #./ (omega_gamma_value ./ γ_value)
            jvalue_mean_rotation_angle = mean(entry.mean_rotation_angles for entry in matching_omega_gamma_data) #./ (omega_gamma_value ./ γ_value)
            jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)

            
            # Append values using push!
            push!(loop_mean_aspect_ratio_list, jvalue_mean_aspect_ratio)
            push!(loop_mean_rotation_angles, jvalue_mean_rotation_angle)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        loop_mean_aspect_ratio_list = loop_mean_aspect_ratio_list[valid_indices]
        loop_mean_rotation_angles = loop_mean_rotation_angles[valid_indices]
        loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]

        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$\\hat{P} = %.4f, \\hat{\\gamma} = %.2f\$", pressure_value, plot_gamma)

        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        mean_aspect_ratio = $(loop_mean_aspect_ratio_list);
        mean_rotation_angles = $(loop_mean_rotation_angles);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);

        % Plot Attenuation
        loglog(ax_attenuation, omega_gamma, mean_attenuation_x, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
        % Plot Rotation Angle
        plot(ax_rotation, omega_gamma, mean_rotation_angles, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        
        % Plot Aspect Ratio
        plot(ax_aspect_ratio, omega_gamma, mean_aspect_ratio, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        """
    end

    # Add legends to the plots
    mat"""
    legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    %legend(ax_rotation, 'show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    %legend(ax_aspect_ratio, 'show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end
#-----------------------------------------------------------------------------------------------------------

function plotGausPressDep(simulation_data, omega_value; plot=true)
    loop_mean_attenuation_list = []
    matching_omega_data = filterDataGaus(simulation_data, omega_value, :omega)

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_omega_data])) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    if plot
        mat"""
        figure_attenuation = figure;
        xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log')
        grid on;
        box on;
        """
    end
        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];
        loop_std_attenuation_list = Float64[]
    for pressure_value in pressure_list

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_omega_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression



        loop_mean_alphaoveromega = mean(filter(x -> x > 0, [entry.attenuation for entry in matching_pressure_data])) ./ omega_value
        loop_std_alphaoveromega =  std(entry.attenuation for entry in matching_pressure_data) ./ omega_value
        push!(loop_std_attenuation_list, loop_std_alphaoveromega)

        push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
    end
    
    if plot
        # Transfer data to MATLAB
        mat"""
        x = $(pressure_list);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        figure(figure_attenuation);
        set(gca, 'Yscale', 'log');
        loglog(x, mean_attenuation_x, '-o')
        hold on
        loglog([x(2), x(end)], [mean_attenuation_x(2), mean_attenuation_x(2)*(x(end)/x(2))^-.75], '--')
        grid on
        box on
        """
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return loop_mean_attenuation_list, pressure_list 
end

function plotAmpRatio(simulation_data, γ_value)
   
    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
    plot_gamma = γ_value
    gamma_value = γ_value

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry

    # Limit range to data
    upper_limit_line_x = [1*γ_value; 1*γ_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*γ_value; .1*γ_value]
    lower_limit_line_y = [1E-5; 1]

    # # Start MATLAB session
    # mat"""
    # figure_main = figure;
    # %tiled_main = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'none'); % 2 rows, 1 column

    # % Axes for Attenuation
    # ax_attenuation = nexttile;
    # hold(ax_attenuation, 'on');
    # % xlabel(ax_attenuation, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    # ylabel(ax_attenuation, '\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}}\$', "FontSize", 20, "Interpreter", "latex");
    # set(ax_attenuation, 'XScale', 'log');
    # set(ax_attenuation, 'YScale', 'log')
    # set(get(ax_attenuation, 'ylabel'), 'rotation', 0);
    # grid(ax_attenuation, 'on');
    # box(ax_attenuation, 'on');
    # %plot(ax_attenuation, $(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$');
    # %plot(ax_attenuation, $(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$');
    # set(ax_attenuation, 'XTickLabel', []);

    # % Axes for Energy
    # ax_energy = nexttile;
    # hold(ax_energy, 'on');
    # xlabel(ax_energy, '\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    # ylabel(ax_energy, '\$ \\overline{\\frac{A_{\\parallel}}{A_{\\perp}}} \$', "FontSize", 15, "Interpreter", "latex");
    # set(ax_energy, 'XScale', 'log');
    # set(ax_energy, 'YScale', 'log')
    # set(get(ax_energy, 'ylabel'), 'rotation', 0);
    # grid(ax_energy, 'on');
    # box(ax_energy, 'on');
    # %set(ax_energy, 'XTickLabel', []);
    # """

    mat"""
    ax_energy = figure;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$  \\overline{\\frac{A_{\\perp}}{A_{\\parallel}}} \$', "FontSize", 20, "Interpreter", "latex");
    set(get(gca, 'ylabel'), 'rotation', 0);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;
    """

    # get a range for plotting color from 0 to 1
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_E_list = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
        
        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current omega_gamma value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds       
            jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
            E_ratio_list = Float64[]
            seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))
            for k_seed in seed_list
                k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                k_seed_omega = k_seed_data[1].omega
                phase_vector_y = k_seed_data[1].unwrapped_phase_vector_y
                # Wrap the phase vector around 2π
                wrapped_phase = mod.(phase_vector_y, 2π)
                distance_from_wall = k_seed_data[1].initial_distance_from_oscillation_output_y_fft
                # mean_distance = meanDistNeighbor(distance_from_wall, wrapped_phase)
                mean_distance = plotAmp(k_seed_data; plot=false)
                push!(E_ratio_list, mean_distance)
            end
            j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
            push!(loop_mean_E_list, j_E_ratio)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        # valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        # matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        # loop_mean_E_list = loop_mean_E_list[valid_indices]
        # loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]
        @bp
        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$ %.4f\$", pressure_value)

        # Transfer data to MATLAB
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        loop_mean_E_list = $(loop_mean_E_list);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);
        
        % Plot Aspect Ratio
        loglog( omega_gamma, loop_mean_E_list, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        ylim([.01, .6])
        xlim([min(oemga_gamma), 1])
        """
    end

    # Add legends to the plots
    mat"""
    leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    title(leg, "\$ \\hat{P} \$")
    """ 
        
end

function plotPhaseRatio(simulation_data, γ_value) 

    # filter the data based on those that are close to gamma_value
    closest_γ_index = argmin(abs.([idx.gamma for idx in simulation_data] .- γ_value))
    closest_γ_value = simulation_data[closest_γ_index].gamma
    matching_γ_data = filter(entry -> entry.gamma == closest_γ_value, simulation_data)
    plot_gamma = γ_value
    gamma_value = γ_value

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_γ_data])) # goes through each entry of simulation_data and get the P value at that entry

    # Limit range to data
    upper_limit_line_x = [1*γ_value; 1*γ_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*γ_value; .1*γ_value]
    lower_limit_line_y = [1E-5; 1]

    mat"""
    ax_energy = figure;
    %xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ 1-\\cos \\overline{\\sigma}_{\\Delta \\phi_{\\perp}} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;
    """

    # get a range for plotting color from 0 to 1
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each pressure
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(element -> element == pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_γ_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_E_list = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))
        
        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current omega_gamma value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds       
            jvalue_mean_alphaoveromega = mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)
            E_ratio_list = Float64[]
            seed_list = sort(unique([entry.seed for entry in matching_omega_gamma_data]))
            for k_seed in seed_list
                k_seed_data = FilterData(matching_omega_gamma_data, k_seed, :seed)
                k_seed_omega = k_seed_data[1].omega
                phase_vector_y = k_seed_data[1].unwrapped_phase_vector_y
                # Wrap the phase vector around 2π
                wrapped_phase = mod.(phase_vector_y, 2π)
                distance_from_wall = k_seed_data[1].initial_distance_from_oscillation_output_y_fft
                # mean_distance = meanDistNeighbor(distance_from_wall, wrapped_phase)
                if isempty(k_seed_data[1].unwrapped_phase_vector_y)
                    println("Empty y-phase vector for: Pressure $(pressure_value) OmegaGamma $(omega_gamma_value) seed $(k_seed)")
                    continue
                end
                mean_distance = plotPhase(k_seed_data; plot=false)
                mean_distance = isinf(mean_distance) ? NaN : mean_distance # saftey for infinite values
                push!(E_ratio_list, 1-cos(mean_distance))
                # push!(E_ratio_list, mean_distance)
            end
            j_E_ratio = mean(E_ratio_list) # mean of the seeds for a single simulation
            push!(loop_mean_E_list, j_E_ratio)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end


        # This is needed because MATLAB.jl has a hard time escaping \'s
        pressure_label = @sprintf("\$ %.4f\$", pressure_value)

        gamma_val = γ_value
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        loop_mean_E_list = $(loop_mean_E_list);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = $(pressure_label);

        %plot( omega_gamma/$(gamma_val), loop_mean_E_list, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        plot( omega_gamma, loop_mean_E_list, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        """
    end

    # Add legends to the plots
    mat"""
    % legend(ax_attenuation, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex');
    leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    title(leg, "\$ \\hat{P} \$")
    """
end
    

function plotGausWavespeed2d(simulation_data; plot=true)  
    loop_mean_wavespeed_list = []
    pressure_out = Float64[]
    wavespeed_out = Float64[]

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        % loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        loglog($(theory_x), $(theory_x).^.5, 'k:', 'DisplayName', 'Slope 1/2');
        xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\hat{c} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each omega value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, simulation_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_wavespeed_list = Float64[];
        loop_std_wavespeed_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))

        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            loop_mean_wavespeed = mean(filter(x -> x > 0, [entry.wavespeed for entry in matching_omega_data]))
            loop_std_wavespeed =  std(entry.wavespeed for entry in matching_omega_data) 
            # loop_mean_wavespeed = mean(filter(x -> x > 0, [entry.wavespeed for entry in matching_omega_data])) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            # loop_std_wavespeed =  std(entry.wavespeed for entry in matching_omega_data) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            push!(loop_std_wavespeed_list, loop_std_wavespeed)

            # Append values
            push!(loop_mean_wavespeed_list, loop_mean_wavespeed)
        end
        mean_wavespeed = mean(loop_mean_wavespeed_list)
        std_wavespeed = std(loop_mean_wavespeed_list)
        lower_bound = mean_wavespeed - 2 * std_wavespeed
        upper_bound = mean_wavespeed + 2 * std_wavespeed
        cleaned_wavespeed = filter(x -> x >= lower_bound && x <= upper_bound, loop_mean_wavespeed_list)
        mean_cleaned_wavespeed = mean(cleaned_wavespeed)
        println("$pressure_value")
        pressure_out = push!(pressure_out, pressure_value)
        wavespeed_out = push!(wavespeed_out, mean_cleaned_wavespeed)

        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.4f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            x = $(matching_omega_list);
            y = $(loop_mean_wavespeed_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(x, y, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            % errorbar(x, y, $(loop_std_wavespeed_list), '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return pressure_out, wavespeed_out
end

#### Gaus Plots

function plotGausWavenumber2d(simulation_data; plot=true)
    loop_mean_wavenumber_list = []
    wavenumber_out = Float64[]
    pressure_out, wavespeed_out = plotGausWavespeed2d(simulation_data; plot=false) 

    # Get a list of unique input pressures
    # pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    pressure_list = sort(unique(pressure for pressure in pressure_out))
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        % loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        loglog($(theory_x), $(theory_x).^.5, 'k:', 'DisplayName', 'Slope 1/2');
        xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\hat{k} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each omega value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]
        wavespeed_at_pressure = wavespeed_out[idx]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, simulation_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_wavenumber_list = Float64[];
        loop_std_wavenumber_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))

        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            # loop_mean_wavespeed = mean(filter(x -> x > 0, [entry.wavespeed for entry in matching_omega_data])) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter)
            loop_std_wavespeed =  std(entry.wavespeed for entry in matching_omega_data) ./ ( sqrt.(matching_omega_data[1].spring_constant ./ matching_omega_data[1].mass) .* matching_omega_data[1].mean_diameter )
            loop_mean_wavespeed = wavespeed_at_pressure

            loop_mean_wavenumber = omega_value ./ loop_mean_wavespeed 
            loop_std_wavenumber= omega_value ./ loop_std_wavespeed

            push!(loop_std_wavenumber_list, loop_std_wavenumber)

            # Append values
            push!(loop_mean_wavenumber_list, loop_mean_wavenumber)
        end


        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.4f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            x = $(matching_omega_list);
            y = $(loop_mean_wavenumber_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(x, y, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            % errorbar(x, y, $(loop_std_wavenumber_list), '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return wavenumber_out, pressure_out
end

function plotGausAttenuation2d(simulation_data; plot=true)  
    loop_mean_attenuation_list = []


    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        % loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        % loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        % loglog($(theory_x), $(theory_x).^.5 ./ 15, 'k:', 'DisplayName', 'Slope 1/2');
        loglog($(theory_x), $(theory_x).^1.5 ./ 15, 'k', 'DisplayName', 'Slope 3/2');
        xlabel('\$\\hat{\\omega}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each omega value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, simulation_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];
        loop_std_attenuation_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))

        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean(filter(x -> x > 0, [entry.attenuation for entry in matching_omega_data])) ./ omega_value
            loop_std_alphaoveromega =  std(entry.attenuation for entry in matching_omega_data) ./ omega_value
            push!(loop_std_attenuation_list, loop_std_alphaoveromega)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
        end


        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.3f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(matching_omega_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(omega_gamma, mean_attenuation_x, '-o', 'Color', marker_color, 'DisplayName', legend_label)
            % errorbar(omega_gamma, mean_attenuation_x, $(loop_std_attenuation_list), '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return loop_mean_attenuation_list
end

function plotGausAttenuationK2d(simulation_data; plot=true)  
    loop_mean_attenuation_list = []
    pressure_out, wavespeed_out = plotGausWavespeed2d(simulation_data; plot=false)

    # Get a list of unique input pressures
    # pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    pressure_list = sort(unique(pressure for pressure in pressure_out))
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        % loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        loglog($(theory_x), $(theory_x).^.5 ./ 15, 'k:', 'DisplayName', 'Slope 1/2');
        loglog($(theory_x), $(theory_x).^1.5 ./ 15, 'k', 'DisplayName', 'Slope 3/2');
        xlabel('\$\\hat{k}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each omega value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, simulation_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];
        loop_std_attenuation_list = Float64[]
        wavenumber_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        # matching_wavenumber_list = sort(unique([entry.wavenumber for entry in matching_pressure_data]))
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))
        # matching_wavenumber_list = wavenumber_out
        
        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            # matching_wavenumber_data = filter(entry -> entry.wavenumber == wavenumber_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            # omega_value = matching_wavenumber_data[1].omega
            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean(filter(x -> x > 0, [entry.attenuation for entry in matching_omega_data]))  ./ omega_value
            loop_std_alphaoveromega =  std(entry.attenuation for entry in matching_omega_data) ./ omega_value
            push!(loop_std_attenuation_list, loop_std_alphaoveromega)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
            loop_wavenumber_ = omega_value / wavespeed_out[idx]
            push!(wavenumber_list, loop_wavenumber_)
        end
        
        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.4f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(wavenumber_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(omega_gamma, mean_attenuation_x, '-|', 'Color', marker_color, 'DisplayName', legend_label)
            % errorbar(omega_gamma, mean_attenuation_x, $(loop_std_attenuation_list), '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return loop_mean_attenuation_list
end

#### Regular Plots

function plot_ωγ_attenuation_2d(simulation_data, gamma_value, mean_diameter; plot=true)
    # Initialize outputs
    matching_omega_gamma_list = []
    loop_mean_attenuation_list = []

    # filter the data based on those that are close to gamma_value
    closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
    closest_gamma_value = simulation_data[closest_gamma_index].gamma
    matching_gamma_data = filter(entry -> entry.gamma == closest_gamma_value, simulation_data)
    plot_gamma = gamma_value

    # Get a list of unique input pressures
    pressure_list = sort(unique([entry.pressure for entry in matching_gamma_data])) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        hold on;
        xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$ \\frac{\\hat{\\alpha}}{\\hat{\\omega}} \$', "FontSize", 20, "Interpreter", "latex");
        set(gca, 'XScale', 'log');
        set(get(gca, 'ylabel'), 'rotation', 0);
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each gamma value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))

        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current pressure value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean_diameter .* mean(entry.alphaoveromega_x for entry in matching_omega_gamma_data)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        loop_mean_attenuation_list = loop_mean_attenuation_list[valid_indices]

        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.3f, \\hat{\\gamma} = %.3f\$", pressure_value, gamma_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(matching_omega_gamma_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            plot_gamma = $(plot_gamma);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(omega_gamma, mean_attenuation_x, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return matching_omega_gamma_list, loop_mean_attenuation_list
end

function plot_ωγ_wavespeed_2d(simulation_data, gamma_value) # Need to fix lgend

    # filter the data based on those that are close to gamma_value
    closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
    closest_gamma_value = simulation_data[closest_gamma_index].gamma
    matching_gamma_data = filter(entry -> entry.gamma == closest_gamma_value, simulation_data)
    plot_gamma = gamma_value

    # Get a list of unique input pressures
    pressure_list = unique([entry.pressure for entry in matching_gamma_data]) # goes through each entry of simulation_data and get the P value at that entry
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = 1 ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(0.5) ./ (1 .+ theory_x.^2);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    # Intialized the plots to loop over
    mat"""
    figure_wavespeed = figure;
    loglog($(theory_x), 1./$(theory_y), 'k', 'DisplayName', '1-D Theory'), hold on
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$\\hat{c} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    """

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each gamma value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_wavespeed_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        matching_omega_gamma_list = sort(unique([entry.omega_gamma for entry in matching_pressure_data]))

        for omega_gamma_value in matching_omega_gamma_list

            # Only look at data for current pressure value
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            @bp
            # Get the mean over all seeds
            loop_mean_wavespeed = mean(entry.:wavespeed_x for entry in matching_omega_gamma_data)
            @bp
            # Append values
            push!(loop_mean_wavespeed_list, loop_mean_wavespeed)
        end

        # Filter data to include only points where omega_gamma <= gamma_value
        valid_indices = matching_omega_gamma_list .<= gamma_value.*2
        matching_omega_gamma_list = matching_omega_gamma_list[valid_indices]
        loop_mean_wavespeed_list = loop_mean_wavespeed_list[valid_indices]
        
        mat"""
        omega_gamma = $(matching_omega_gamma_list);
        mean_wavespeed_x = $(loop_mean_wavespeed_list);
        iloop_pressure_value = $(pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = sprintf('Wavespeed X = %.2f, Gamma = %.2f (Aspect Ratio)', $(pressure_value), $(plot_gamma));
        % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
        pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Wavespeed X)', $(pressure_value), $(plot_gamma));
        
        figure(figure_wavespeed);
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_wavespeed_x, '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
        """
    end

    # Add legends to the plots
    mat"""
    % Add legends to the MATLAB plots
    figure(figure_wavespeed);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end

function plotGuess(simulation_data; plot=true)  # This is expirimental
    loop_mean_attenuation_list = []
    pressure_out, wavespeed_out = plotGausWavespeed2d(simulation_data; plot=false)

    # Get a list of unique input pressures
    # pressure_list = sort(unique([entry.pressure for entry in simulation_data])) # goes through each entry of simulation_data and get the P value at that entry
    pressure_list = sort(unique(pressure for pressure in pressure_out))
    plot_pressure = pressure_list

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);
    # upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    # upper_limit_line_y = [1E-5; 1]
    # lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    # lower_limit_line_y = [1E-5; 1]
    # % plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    # % plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')

    if plot
        # Intialized the plots to loop over
        mat"""
        figure_attenuation = figure;
        % loglog($(theory_x), $(theory_y), 'k', 'DisplayName', '1-D Theory');
        %loglog($(theory_x), $(theory_x) ./ 15, 'k--', 'DisplayName', 'Slope 1');
        hold on;
        %loglog($(theory_x), $(theory_x).^.5 ./ 15, 'k:', 'DisplayName', 'Slope 1/2');
        %loglog($(theory_x), $(theory_x).^1.5 ./ 15, 'k', 'DisplayName', 'Slope 3/2');
        xlabel('\$\\hat{k}\\hat{P}\$', "FontSize", 20, "Interpreter", "latex");
        ylabel('\$\\hat{\\alpha} \$', "FontSize", 20, "Interpreter", "latex");
        set(get(gca, 'ylabel'), 'rotation', 0);
        set(gca, 'XScale', 'log');
        grid on;
        box on;
        """
    end

    # Normalize the gamma values
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    # Create a line for each omega value across all pressure_list
    for pressure_value in pressure_list

        # Assign a color
        idx = findfirst(idx -> idx ==pressure_value, pressure_list) # find the first index that matches
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # Only look at data for current pressure value
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, simulation_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];
        loop_std_attenuation_list = Float64[]
        wavenumber_list = Float64[]

        # Look at a single omega gamma value since each one spans all seeds
        # matching_wavenumber_list = sort(unique([entry.wavenumber for entry in matching_pressure_data]))
        matching_omega_list = sort(unique([entry.omega for entry in matching_pressure_data]))
        # matching_wavenumber_list = wavenumber_out
        
        for omega_value in matching_omega_list

            # Only look at data for current pressure value
            # matching_wavenumber_data = filter(entry -> entry.wavenumber == wavenumber_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            matching_omega_data = filter(entry -> entry.omega == omega_value, matching_pressure_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
            # omega_value = matching_wavenumber_data[1].omega
            # Get the mean over all seeds
            loop_mean_alphaoveromega = mean(filter(x -> x > 0, [entry.attenuation for entry in matching_omega_data]))  #./ omega_value
            loop_std_alphaoveromega =  std(entry.attenuation for entry in matching_omega_data) #./ omega_value
            push!(loop_std_attenuation_list, loop_std_alphaoveromega)

            # Append values
            push!(loop_mean_attenuation_list, loop_mean_alphaoveromega)
            loop_wavenumber_ = omega_value / wavespeed_out[idx] #* pressure_value
            push!(wavenumber_list, loop_wavenumber_)
        end
        
        # This is needed because MATLAB.jl has a hard time escaping \'s
        legend_label = @sprintf("\$ \\hat{P} = %.4f\$", pressure_value)

        if plot
            # Transfer data to MATLAB
            mat"""
            omega_gamma = $(wavenumber_list);
            mean_attenuation_x = $(loop_mean_attenuation_list);
            iloop_pressure_value = $(pressure_value);
            marker_color= $(marker_color);
            legend_label = $(legend_label);
            figure(figure_attenuation);
            set(gca, 'Yscale', 'log');
            plot(omega_gamma, mean_attenuation_x, '-|', 'Color', marker_color, 'DisplayName', legend_label)
            xline($(pressure_value) * .5)
            yline(.12)
            % errorbar(omega_gamma, mean_attenuation_x, $(loop_std_attenuation_list), '-o','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', legend_label);
            """
        end
    end
    if plot
        # Add legends to the plots
        mat"""
        % Add legends to the MATLAB plots
        figure(figure_attenuation);
        legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
        """
    end
    return loop_mean_attenuation_list
end

# Single Simulaiton plots

function plotAmp(filtered_data; plot=true, shear=false)
    x_parra = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    y = filtered_data[1].amplitude_vector_x
    coeffs = fitLogLine(x_parra,y)
    yIntercept_amp_x = coeffs[1]
    slope_amp_x = coeffs[2]
    y_x = y
    para_color = 'b'

    if shear==true
       x_parra = filtered_data[1].initial_distance_from_oscillation_output_y_fft
       y = filtered_data[1].amplitude_vector_y
       coeffs = fitLogLine(x_parra,y)
       yIntercept_amp_x = coeffs[1]
       slope_amp_x = coeffs[2]
       y_x = y
       para_color = 'r'
    end
    if plot==true
        display_name = @sprintf("\$ A_{||}(x) \$")
        mat"""
        figure
        scatter($(x_parra), $(y), 'Color', $(para_color), "DisplayName", $(display_name))
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$ A(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        """
    end

    x_perp = filtered_data[1].initial_distance_from_oscillation_output_y_fft
    y = filtered_data[1].amplitude_vector_y
    coeffs = fitLogLine(x_perp,y)
    yIntercept_amp_y = coeffs[1]
    slope_amp_y = coeffs[2]
    perp_color = 'r'
    if shear == true
        x_perp = filtered_data[1].initial_distance_from_oscillation_output_x_fft
        y = filtered_data[1].amplitude_vector_x
        coeffs = fitLogLine(x_perp,y)
        yIntercept_amp_y = coeffs[1]
        slope_amp_y = coeffs[2]
        perp_color = 'b'
    end
    if plot == true
        display_name = @sprintf("\$ A_\\perp(x) \$")
        mat"""
        scatter($(x_perp), $(y), 'Color', $(perp_color), "DisplayName", $(display_name))
        plot($(x_perp), exp($(yIntercept_amp_y) + $(slope_amp_y) .* $(x_perp)), 'HandleVisibility', 'off')
        plot($(x_parra), exp($(yIntercept_amp_x) + $(slope_amp_x) .* $(x_parra)), 'HandleVisibility', 'off')
        set(gca, 'YScale', 'log')
        grid on
        legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
        """
    end
    
    return exp(yIntercept_amp_y) / exp(yIntercept_amp_x) # Can't do  mean(y ./ y_x), because not same length
end

function plotPhase(filtered_data; plot=true, shear=false)
    distance_y = filtered_data[1].initial_distance_from_oscillation_output_y_fft
    distance_x = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    phase_y = filtered_data[1].unwrapped_phase_vector_y
    phase_x = filtered_data[1].unwrapped_phase_vector_x
    phase_y = mod.(phase_y, 2π)
    phase_x = mod.(phase_x, 2π)
    scatter_x = meanPhaseDev(distance_x, phase_x, 1)
    scatter_y = meanPhaseDev(distance_y, phase_y, 1)
    para_color = "b"
    perp_color = "r"

    if shear==true
        # switch the vaalues of scatter_x and scatter_y
        temp = scatter_x
        scatter_x = scatter_y
        scatter_y = temp
        temp2 = distance_x
        distance_x = distance_y
        distance_y = temp2
        temp3 = phase_x
        phase_x = phase_y
        phase_y = temp3
        para_color = "r"
        perp_color = "b"
    end

    if plot==true
        mat"""
        figure
        scatter($(distance_x), $(phase_x), 'Color', $(para_color), "DisplayName", "\$ \\phi_{||} \$")
        hold on
        scatter($(distance_y), $(phase_y), 'Color', $(para_color), "DisplayName", "\$ \\phi_{\\perp} \$")
        grid on
        box on
        set(gca,'YTick', [0, pi, 2*pi], 'YTickLabel', {'0', ' \$ \\pi \$', '\$ 2\\pi \$'}, 'TickLabelInterpreter', 'latex');
        ylabel("\$ \\Delta \\phi \$", "Interpreter", 'latex', "FontSize", 15)
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        legend('Interpreter', 'latex')
        ylim([0,2*pi])
        """
    end

    return  scatter_y
end

function combinePlots()
    mat"""
        f1 = 'fig1.fig';
        f2 = 'fig2.fig';

        % Open the figures
        fig1 = openfig(f1);
        fig2 = openfig(f2);

        % get axes handles - this assumes there is only 1 axes per figure!
        fig1ax = gca(fig1);
        fig2ax = gca(fig2);
        leg1 = findobj(fig1,'Type','legend');
        leg2 = findobj(fig2,'Type','legend');

        % Get axis children
        fig1axChildren = get(fig1ax,'Children');
        fig2axChildren = get(fig2ax,'Children');

        % Create new fig, copy items from fig 1
        % This will maintain set properties such as color
        figFinal = figure();
        ax = axes(figFinal);
        h1 = copyobj(fig1axChildren, ax);


        % Copy items from fig 2
        h2 = copyobj(fig2axChildren, ax);
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        grid on

        % Add legend to same location as the legend in fig2 
        % but only include objects with a defined DisplayName
        h = [h2;h1];
        hasDisplayName = ~cellfun('isempty',get(h,'DisplayName'));
        legend(ax, h(hasDisplayName),'Location', leg1.Location, 'Interpreter', 'Latex')

        % Copy axis labels
        xlabel(ax, fig2ax.XLabel.String, 'Interpreter', 'Latex', 'FontSize', 20)
        ylabel(ax, fig2ax.YLabel.String, 'Interpreter', 'Latex', 'FontSize', 20)
        set(get(gca, 'ylabel'), 'rotation', 0);
    """
end

function plotAmpTiled(filtered_data_a; plot=true)
    x_parra_a = filtered_data_a[1].initial_distance_from_oscillation_output_x_fft
    y_a = filtered_data_a[1].amplitude_vector_x

    if plot==true
        display_name = @sprintf("\$ A_{||}(x) \$")
        mat"""
        figure_main = figure;
        tiled_main = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'none'); % 2rows, 1 column
        figure_amp = nexttile
        scatter(figure_amp, $(x_parra_a), $(y_a), "DisplayName", $(display_name))
        hold(figure_amp, 'on')
        set(figure_amp, 'YScale', 'log')
        grid(figure_amp, 'on')
        ylabel(figure_amp, "\$ A(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(figure_amp, 'ylabel'), 'rotation', 0);
        set(figure_amp, 'XTickLabel', []);
        box on
        hold on 
        """
    end

    x_perp_a = filtered_data_a[1].initial_distance_from_oscillation_output_y_fft
    y_a = filtered_data_a[1].amplitude_vector_y
    if plot == true
        display_name = @sprintf("\$ A_\\perp(x) \$")
        mat"""
        scatter(figure_amp, $(x_perp_a), $(y_a), "DisplayName", $(display_name))
        set(figure_amp, gca, 'YScale', 'log')
        
        grid on
        """
    end
    filtered_data = filtered_data_a # i'm lazy
    distance_y = filtered_data[1].initial_distance_from_oscillation_output_y_fft
    distance_x = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    phase_y = filtered_data[1].unwrapped_phase_vector_y
    phase_x = filtered_data[1].unwrapped_phase_vector_x
    phase_y = mod.(phase_y, 2π)
    phase_x = mod.(phase_x, 2π)
    # scatter_x = meanDistNeighbor(distance_x, phase_x)
    # scatter_y = meanDistNeighbor(distance_y, phase_y)
    scatter_x = meanDeltaYNeighbor(distance_x, phase_x)
    scatter_y = meanDeltaYNeighbor(distance_y, phase_y)
    if plot==true
        mat"""
        figure_phase = nexttile
        scatter(figure_phase, $(distance_x), $(phase_x), "DisplayName", "\$ \\hat{x} \$")
        hold on
        scatter(figure_phase, $(distance_y), $(phase_y), "DisplayName", "\$ \\hat{y} \$")
        grid(figure_phase, 'on')
        box(figure_phase, 'on')
        set(figure_phase,'YTick', [0, pi, 2*pi], 'YTickLabel', {'0', ' \$ \\pi \$', '\$ 2\\pi \$'}, 'TickLabelInterpreter', 'latex');
        ylabel(figure_phase, "\$ \\Delta \\phi \$", "Interpreter", 'latex', "FontSize", 15)
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(figure_phase, 'ylabel'), 'rotation', 0);
        legend(figure_phase, 'Interpreter', 'latex', "FontSize", 15)
        ylim(figure_phase, [0,2*pi])
        """
    end

# This is an artifact, but dont' want to delete yet because might actuall use 
    # x_parra_b = filtered_data_b[1].initial_distance_from_oscillation_output_x_fft
    # y_b = filtered_data_b[1].amplitude_vector_x

    # if plot==true
    #     display_name = @sprintf("\$ A_{||}(x) \$")
    #     mat"""
    #     figure
    #     scatter($(x_parra_b), $(y_b), "DisplayName", $(display_name))
    #     hold on
    #     set(gca, 'YScale', 'log')
    #     grid on
    #     xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
    #     ylabel("\$ A(x) \$", "Interpreter", 'latex', "FontSize", 15)
    #     set(get(gca, 'ylabel'), 'rotation', 0);
    #     box on
    #     hold on 
    #     """
    # end

    # x_perp_b = filtered_data_b[1].initial_distance_from_oscillation_output_y_fft
    # y_b = filtered_data_b[1].amplitude_vector_y
    # if plot == true
    #     display_name = @sprintf("\$ A_\\perp(x) \$")
    #     mat"""
    #     scatter($(x_perp_b), $(y_b), "DisplayName", $(display_name))
    #     set(gca, 'YScale', 'log')
    #     grid on
    #     legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);
    #     """
    # end

end