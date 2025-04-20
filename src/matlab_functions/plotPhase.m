function scatter_y = plotPhase(data, options)

    arguments
        data (1,1) struct % Required: scalar struct with simulation data
        options.plot (1,1) logical = true % Optional: flag to generate plots
        options.shear (1,1) logical = false % Optional: flag for shear vs compression for parallel component
        options.z (1,1) logical = false % plots z direction instead
    end

    if options.z
        distance_y = data.initial_distance_from_oscillation_output_z_fft{1};
        distance_x = data.initial_distance_from_oscillation_output_x_fft{1};
        phase_y = data.unwrapped_phase_vector_z{1};
        phase_x = data.unwrapped_phase_vector_x{1};
        phase_y = mod(phase_y, 2*pi);
        phase_x = mod(phase_x, 2*pi);
        scatter_x = meanPhaseDev(distance_x, phase_x,1);
        scatter_y = meanPhaseDev(distance_y, phase_y, 1);
    else
        distance_y = data.initial_distance_from_oscillation_output_y_fft{1};
        distance_x = data.initial_distance_from_oscillation_output_x_fft{1};
        phase_y = data.unwrapped_phase_vector_y{1};
        phase_x = data.unwrapped_phase_vector_x{1};
        phase_y = mod(phase_y, 2*pi);
        phase_x = mod(phase_x, 2*pi);
        scatter_x = meanPhaseDev(distance_x, phase_x,1);
        scatter_y = meanPhaseDev(distance_y, phase_y, 1);
    end

    if options.shear
        % # switch the vaalues of scatter_x and scatter_y. There's probably some O(n) leet code answer for this
        temp = scatter_x
        scatter_x = scatter_y
        scatter_y = temp
        temp2 = distance_x
        distance_x = distance_y
        distance_y = temp2
        temp3 = phase_x
        phase_x = phase_y
        phase_y = temp3
    end

    if options.plot
        figure
        scatter(distance_x, phase_x, "DisplayName", "$ \phi_{||} $")
        hold on
        scatter(distance_y, phase_y, "DisplayName", "$ \phi_{\perp} $")
        grid on
        box on
        set(gca,'YTick', [0, pi, 2*pi], 'YTickLabel', {'0', ' $ \pi $', '$ 2\pi $'}, 'TickLabelInterpreter', 'latex');
        ylabel("$ \Delta \phi $", "Interpreter", 'latex', "FontSize", 15)
        xlabel("$ x $", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        legend('Interpreter', 'latex')
        ylim([0,2*pi])
    end
end