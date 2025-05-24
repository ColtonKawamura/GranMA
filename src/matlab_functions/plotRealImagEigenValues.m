function plotRealImagEigenValues(plotData, pressure_list, damping_list, options)
% plotRealImagEigenValues - Plot real and imaginary parts of eigenvalues
%
% Syntax: plotRealImagEigenValues(plotData, pressure_list, damping_list, options)
%
% Inputs:
%    plotData - data structure containing eigenvalues and other properties
%    pressure_list - array of pressures to plot
%    damping_list - array of damping constants to plot
%    options - structure containing optional parameters
%        .loglog - whether to use log scale for the axes (default: false)
%        .semilogx - whether to use semi-log scale for x-axis (default: false)
%        .semilogy - whether to use semi-log scale for y-axis (default: false)
%        .scaling - whether to scale the axes (default: false)
%
% Outputs:
%    None
%
% Example:
%    plotRealImagEigenValues(plotData, [0.1, 0.2], [0.5, 1.0], struct('loglog', true))
%    This will plot the real and imaginary parts of eigenvalues for the specified pressures and damping constants
%    with log-log scaling.
%    and the specified options.
    arguments
        plotData % data structure containing eigenvalues and other properties
        pressure_list (1,:) double % list of pressures to plot
        damping_list (1,:) double % list of damping constants to plot
        options.loglog (1,1) logical = false % whether to use log scale for the axes
        options.semilogx (1,1) logical = false % whether to use semi-log scale for x-axis
        options.semilogy (1,1) logical = false % whether to use semi-log scale for y-axis
        options.scaling (1,1) logical = false % whether to scale the axes
    end

    figure
    for i = 1:length(pressure_list)
        pressure = pressure_list(i);

        for j = 1:length(damping_list)
            damping_constant = damping_list(j);
            dataPressureDamping = filterData(plotData, 'pressure', pressure,  'damping', damping_constant);
            realEigenValues = real(dataPressureDamping.eigenValues{1});
            imagEigenValues = imag(dataPressureDamping.eigenValues{1});
            keepIdx = imagEigenValues >= 0; % keep the positive eigenvalues
            realEigenValues = realEigenValues(keepIdx); % imag part carries the frequency, abs() because QZ sovler does weird things
            imagEigenValues = imagEigenValues(keepIdx);
            
            Lx = dataPressureDamping.Lx(1);
            Ly = dataPressureDamping.Ly(1);

            if isscalar(pressure_list)
                % [~, marker_color] = normVarColor(damping_list, damping_constant, 1);
                marker_color = [0,0,1];
            else
                [~, marker_color] = normVarColor(pressure_list, pressure, 1);
            end

            if isscalar(damping_list)
                markerSize = 6;
            else
                markerSize = exp(dataPressureDamping.damping/max(damping_list))*3;
            end

            pressureLabel = sprintf('$ %.4f, %.4f $', dataPressureDamping.pressure, dataPressureDamping.damping); 
            pressureValue = dataPressureDamping.pressure;

            if options.scaling
                plot(imagEigenValues./sqrt(pressureValue), -realEigenValues./(sqrt(pressureValue).*damping_constant), 'o', 'MarkerSize', markerSize , 'MarkerEdgeColor', marker_color, 'Color', marker_color, 'DisplayName', pressureLabel);
                ylabel('Re$(\lambda)/(\sqrt{\hat{P}\hat{\gamma}}}$', 'Interpreter', 'latex', 'FontSize', 20)
                xlabel('Im$(\lambda)/\sqrt{\hat{P}}$', 'Interpreter', 'latex', 'FontSize', 20)
            else
                plot(imagEigenValues, -realEigenValues./damping_constant, 'o', 'MarkerSize', markerSize , 'MarkerEdgeColor', marker_color, 'Color', marker_color, 'DisplayName', pressureLabel);
                ylabel('Re$(\lambda)/\hat{\gamma}$', 'Interpreter', 'latex', 'FontSize', 20)
                xlabel('Im$(\lambda)$', 'Interpreter', 'latex', 'FontSize', 20)
            end

            if options.loglog
                set(gca, 'XScale', 'log', 'YScale', 'log');
            elseif options.semilogx
                set(gca, 'XScale', 'log');
            elseif options.semilogy
                set(gca, 'YScale', 'log');
            end
            hold on

            title(sprintf('$L_x$ by $L_y$: %.2f by %.2f', Lx, Ly), 'Interpreter', 'latex', 'FontSize', 16);
            % set(gca, "XScale", "log")
            % set(gca, "YScale", "log")
            grid on; 
            hold on;
        end
        legend('show', 'Interpreter', 'latex');
        leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
        title(leg, "$  \hat{P}, \hat{\gamma} $")
        ax = gca;
        ax.FontSize = 20;
    end
end