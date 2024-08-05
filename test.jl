using GLMakie
using Statistics
#  This does not work in the REPL
x = collect(1:1:4)
y = 2*x

fig = Figure()

ax = Axis(fig[1, 1],
    title = L"\textrm{tst} \alpha",
    xlabel = "Time (seconds)",
    ylabel = L"\hat{\alpha}",
)
scatter!(ax, x, y)
fig

function main()
    x = collect(1:1:4)
    y = 2*x

    fig = Figure()

    ax = Axis(fig[1, 1],
        title = L"\textrm{tst} \alpha",
        xlabel = "Time (seconds)",
        ylabel = L"\hat{\alpha}",
    )
    scatter!(ax, x, y)
    DataInspector(fig)
end

function adv()
    seconds = 0:0.1:2
    measurements = [8.2, 8.4, 6.3, 9.5, 9.1, 10.5, 8.6, 8.2, 10.5, 8.5, 7.2,
        8.8, 9.7, 10.8, 12.5, 11.6, 12.1, 12.1, 15.1, 14.7, 13.1]
        f = Figure()
        ax = Axis(f[1, 1],
            title = "Experimental data and exponential fit",
            xlabel = "Time (seconds)",
            ylabel = "Value",
        )
        scatter!(
            ax,
            seconds,
            measurements,
            color = :tomato,
            label = "Measurements"
        )
        lines!(
            ax,
            seconds,
            exp.(seconds) .+ 7,
            color = :tomato,
            linestyle = :dash,
            label = "f(x) = exp(x) + 7",
        )
        axislegend(position = :rb)
        f

end

function mult()
    test1 = collect(1:1:5)
    test2 = test1.+1
    test3 = test1.+2
    f = Figure()
    ax = Axis(f[1, 1],
        title = L"Experimental data and exponential fit",
        xlabel = "Time (seconds)",
        ylabel = L"\alpha",
        )

    for idx = eachindex(test1)
        testplot = test1.+idx
        lines!(ax, test1, testplot, color=:tomato, label=L"\alpha = %$(idx)")
    end
    axislegend(position = :rb)
    f
end



function testplot_ωγ_attenuation_2d(data_frame, gamma_value)

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

    # Define the plot limits to match the 1D theory plot curves
    theory_x = collect(3E-4:1E-5:3)
    theory_y = theory_x ./ sqrt(2) .* ((1 .+ theory_x.^2) .* (1 .+ sqrt.(1 .+ theory_x.^2))).^(-0.5);

    fig = Figure()
    ax = Axis(fig[1,1], 
        xlabel=L"\hat{\omega}\hat{\gamma}",
        ylabel=L"\frac{\hat{\alpha}}{\hat{\omega}}", 
        xscale = log10,
        yscale = log10)
    lines!(ax,theory_x, theory_y)

    # Normalize the gamma values
    normalized_variable = (plot_pressure .- minimum(plot_pressure)) ./ (maximum(plot_pressure) .- minimum(plot_pressure))

    # Create a line for each gamma value across all pressure_list
    for idx in eachindex(plot_pressure)

        marker_color = RGBf(normalized_variable[idx], 0, 1-normalized_variable[idx])

        # For idx, only show current pressure data
        iloop_pressure_value = plot_pressure[idx]
        iloop_pressure_index = in.(filtered_data_frame.input_pressure, Ref(iloop_pressure_value))
        iloop_combined_index = iloop_pressure_index
        iloop_data_frame = filtered_data_frame[iloop_combined_index, :]

        # Initizalized vectors for just this pressure
        loop_mean_attenuation_list = Float64[];
        loop_omegagamma_plot_list = Float64[];

        # Look at a single omega gamma value since each one spans all seeds
        iloop_omega_gamma_list = sort(unique(iloop_data_frame.omegagamma))

        for jdx in eachindex(iloop_omega_gamma_list)

            # Get the idex for the current omega gamma value
            matching_jdx = in.(iloop_data_frame.omegagamma, Ref(iloop_omega_gamma_list[jdx]))
            jloop_data_frame = iloop_data_frame[matching_jdx,:]

            # get the mean over all seeds
            jvalue_mean_alphaoveromega = mean(jloop_data_frame.alphaoveromega)

        # Append values only if they are positive
        if jvalue_mean_alphaoveromega >= 0
            push!(loop_mean_attenuation_list, jvalue_mean_alphaoveromega)
            push!(loop_omegagamma_plot_list, iloop_omega_gamma_list[jdx])
        end
        end
        scatterlines!(ax, loop_omegagamma_plot_list, loop_mean_attenuation_list, color=marker_color, label=L"\hat{P} = %$(iloop_pressure_value)")
    end
    axislegend(position = :rb)
    fig

end
