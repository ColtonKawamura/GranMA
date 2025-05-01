% filepath: /Users/coltonkawamura/Documents/GranMA/src/matlab_functions/plotAmp.m
function ampRatio = plotAmp(data, options)
    % Plots amplitude decay and calculates the ratio of perpendicular to parallel
    % amplitude intercepts from exponential fits.
    %
    % Args:
    %   filtered_data: A scalar struct containing simulation data fields like
    %                  initial_distance_from_oscillation_output_x_fft,
    %                  amplitude_vector_x, etc.
    %   options: Name-value pairs
    %       plot (logical, default true): Generate plots.
    %       shear (logical, default false): Use shear data (y) instead of compression (x) for parallel component.
    %       transverse_axis (string, default "y"): Which axis ('y' or 'z') to use for the perpendicular component.
    %
    % Returns:
    %   amp_ratio: The ratio exp(yIntercept_amp_perp) / exp(yIntercept_amp_parra)

    arguments
        data (1,1) struct % Required: scalar struct with simulation data
        options.plot (1,1) logical = true % Optional: flag to generate plots
        options.shear (1,1) logical = false % Optional: flag for shear vs compression for parallel component
        options.z (1,1) logical = false % plots z direction instead
    end

    xPara = data.initial_distance_from_oscillation_output_x_fft{1};
    y = data.amplitude_vector_x{1};
    coeffs = fitLogLine(xPara, y);
    yInterecept_xAmp = coeffs(1);
    slopeAmpX = coeffs(2);
    y_x = y;

    if options.plot
        figure
        scatter(xPara, y, "DisplayName", "$A_{||}$")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("$ x $", "Interpreter", 'latex', "FontSize", 20)
        ylabel("$ A(x) $", "Interpreter", 'latex', "FontSize", 20)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
    end

    if options.z
        xPerp = data.initial_distance_from_oscillation_output_z_fft{1};
        y = data.amplitude_vector_z{1};
        coeffs = fitLogLine(xPerp, y);
        yintercept_yAmp = coeffs(1);
        slopeAmpY = coeffs(2);
    else
        xPerp = data.initial_distance_from_oscillation_output_y_fft{1};
        y = data.amplitude_vector_y{1};
        coeffs = fitLogLine(xPerp, y);
        yintercept_yAmp = coeffs(1);
        slopeAmpY = coeffs(2);
    end

    if options.plot
        scatter(xPerp, y, "DisplayName", "$ A_\perp(x) $")
        plot(xPerp, exp(yintercept_yAmp + slopeAmpY .* xPerp), 'HandleVisibility', 'off')
        plot(xPara, exp(yInterecept_xAmp + slopeAmpX .* xPara), 'HandleVisibility', 'off')
        set(gca, 'YScale', 'log')
        grid on
        legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 20);
        title(sprintf("$\\hat{P}$: %g, $\\hat{\\gamma}$: %g, $\\hat{\\omega}$: %g", data.pressure(1), data.gamma(1), data.omega(1)), 'Interpreter', 'latex','FontSize', 20);
    end

    ampRatio = exp(yintercept_yAmp) / exp(yInterecept_xAmp);
end

function coeffs = fitLogLine(x, y)
    % Simple example: fits log(y) = p(1) + p(2)*x
    p = polyfit(x, log(y), 1);
    coeffs = [p(2), p(1)]; % Return as [intercept, slope] consistent with Julia version? Check fitLogLine definition
end