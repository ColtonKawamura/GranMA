module GranMA

using DataFrames
using CSV
using Plots
using LaTeXStrings
using Debugger # REPL: Debugger.@run function(); @bp
using MATLAB
using Statistics

export plot_ellipse_pressure_2d, random_aspect_ratio_check_2d

function main()
    gamma_value = 1
    plot_ellipse_pressure(gamma_value)
    # common_rows = random_aspect_ratio_check()
end

function plot_ellipse_pressure_2d(gamma_value)

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

    # Start MATLAB session
    mat"""
    % Initialize the MATLAB figure for aspect ratio
    figure_aspect_ratio = figure;
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\overline{\\frac{b}{a}} \$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;

    % Initialize the MATLAB figure for rotation angle
    figure_rotation_angle = figure;
    hold on;
    xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\left|\\overline{\\theta} \\right|\$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;
    """

    # Normalize the gamma values
    normalized_variable = (plot_pressure .- minimum(plot_pressure)) ./ (maximum(plot_pressure) .- minimum(plot_pressure))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = RGB(normalized_variable[idx], 0, 1-normalized_variable[idx])
        color_array = [red(marker_color), green(marker_color), blue(marker_color)]

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
        pressure_label = sprintf('Pressure = %.2f, Gamma = %.2f (Aspect Ratio)', $(iloop_pressure_value), $(plot_gamma));
        % pressure_label = "\$ \\alpha x^2 \$ = iloop_pressure_value, \\gamma = plot_gamma \\mathrm{(Attenuation)}"
        pressure_label2 = sprintf('Pressure = %.2f, Gamma = %.2f (Attenuation)', $(iloop_pressure_value), $(plot_gamma));
        color = $(color_array);
        
        figure(figure_aspect_ratio);
        yyaxis left;
        plot(omega_gamma, mean_aspect_ratio, 'o-', 'MarkerFaceColor', color, 'Color', color, 'DisplayName', pressure_label);
        yyaxis right;
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_attenuation_x, '--x','MarkerFaceColor', color, 'Color', color, 'DisplayName', pressure_label2);
        
        figure(figure_rotation_angle);
        yyaxis left;
        plot(omega_gamma, mean_rotation_angles, 'o-', 'MarkerFaceColor', color, 'Color', color, 'DisplayName', pressure_label);
        yyaxis right;
        set(gca, 'Yscale', 'log');
        plot(omega_gamma, mean_attenuation_x, '--x','MarkerFaceColor', color, 'Color', color, 'DisplayName', pressure_label2);
        """
    end

    # Add legends to the plots
    mat"""
    % Add legends to the MATLAB plots
    figure(figure_aspect_ratio);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');

    figure(figure_rotation_angle);
    legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex');
    """
end

function random_aspect_ratio_check_2d()
    # pick a random row 
    random_row = data_frame[rand(1:nrow(data_frame)) ,:]

    # seeds have input_pressure and omega gamma in common
    common_rows = data_frame[in.(data_frame.input_pressure, random_row.input_pressure) .& in.(data_frame.omegagamma, random_row.omegagamma),:]

    return common_rows    
end

end