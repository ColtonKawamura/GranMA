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

function plot_ellipse_ωγ_2d(data_frame, gamma_value, flag::Any) # using MATLAB

    # Define parameters to plot
    pressure_list = sort(unique(data_frame.input_pressure))
    plot_pressure = pressure_list
    gamma_list = sort(unique(data_frame.gamma))
    closest_gamma_match_index = argmin(abs.(gamma_list .- gamma_value))
    plot_gamma = gamma_list[closest_gamma_match_index]

    # Filter the table to only those data
    matching_gamma_index = in.(data_frame.gamma, Ref(plot_gamma))
    combined_index = matching_gamma_index
    filtered_data_frame = data_frame[combined_index, :]

    # Limit range to data
    upper_limit_line_x = [1*gamma_value; 1*gamma_value]
    upper_limit_line_y = [1E-5; 1]
    lower_limit_line_x = [.1*gamma_value; .1*gamma_value]
    lower_limit_line_y = [1E-5; 1]

    # Start MATLAB session
    mat"""
    % Initialize the MATLAB figure for aspect ratio
    figure_aspect_ratio = figure;
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\overline{\\frac{b}{a}} \$', "FontSize", 20, "Interpreter", "latex");
    plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;

    % Initialize the MATLAB figure for rotation angle
    figure_rotation_angle = figure;
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$\\overline{  \\left| \\theta \\right|}\$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    """

    # Normalize the gamma values
    normalized_variable = (log.(plot_pressure) .- minimum(log.(plot_pressure))) ./ (maximum(log.(plot_pressure)) .- minimum(log.(plot_pressure)))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]

        # For idx, only show current pressure data
        iloop_pressure_value = plot_pressure[idx]
        iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
        iloop_combined_index = iloop_pressure_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_aspect_ratio_list = Float64[];
        loop_mean_rotation_angles = Float64[];
        loop_mean_attenuation_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

        for jdx in eachindex(iloop_omega_gamma_list)

            # Get the idex for the current omega gamma value
            matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
            jloop_data_frame = iloop_data_frame[matching_jdx,:]

            # get the mean over all seeds
            jvalue_mean_aspect_ratio = mean(jloop_data_frame.mean_aspect_ratio)
            jvalue_mean_rotation_angle = mean(jloop_data_frame.mean_rotation_angles)
            jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)
            @bp
            # Append values using push!
            push!(loop_mean_aspect_ratio_list, jvalue_mean_aspect_ratio)
            push!(loop_mean_rotation_angles, jvalue_mean_rotation_angle)
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
        end
        
        # Transfer data to MATLAB
        mat"""
        % Transfer data to MATLAB
        omega_gamma = $(iloop_omega_gamma_list);
        mean_aspect_ratio = $(loop_mean_aspect_ratio_list);
        mean_rotation_angles = $(loop_mean_rotation_angles);
        mean_attenuation_x = $(loop_mean_attenuation_list);
        iloop_pressure_value = $(iloop_pressure_value);
        plot_gamma = $(plot_gamma);
        marker_color = $(marker_color);
        pressure_label = sprintf('Pressure = %.2f, Gamma = %.2f (Aspect Ratio)', $(iloop_pressure_value), $(plot_gamma));
        % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
        pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Attenuation)', $(iloop_pressure_value), $(plot_gamma));

        
        figure(figure_aspect_ratio);
        yyaxis left;
        plot(omega_gamma, mean_aspect_ratio, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        yyaxis right;
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_attenuation_x, '--x','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
        
        figure(figure_rotation_angle);
        yyaxis left;
        plot(omega_gamma, mean_rotation_angles, 'o-', 'MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label);
        yyaxis right;
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_attenuation_x, '--x','MarkerFaceColor', marker_color, 'Color', marker_color, 'DisplayName', pressure_label2);
        """
    end

    # Add legends to the plots
    mat"""
    % Add legends to the MATLAB plots
    figure(figure_aspect_ratio);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');

    figure(figure_rotation_angle);
    plot($(upper_limit_line_x), $(upper_limit_line_y), 'k', 'DisplayName', '\$ \\omega_0 \$')
    plot($(lower_limit_line_x), $(lower_limit_line_y), 'b', 'DisplayName', '\$ .1 \\omega_0 \$')
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end