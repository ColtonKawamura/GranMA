using Plots

amplitude_x = 1.0
amplitude_y = 0.5
phase_x = 0.0
phase_y = π/2
omega = 2π
theta = π / 6
t_max = 2.0
n_frames = 100

t_values = range(0, t_max, length=n_frames)

function get_position(t)
    x = amplitude_x * cos(omega * t + phase_x)
    y = amplitude_y * cos(omega * t + phase_y)
    x_rot = cos(theta) * x - sin(theta) * y
    y_rot = sin(theta) * x + cos(theta) * y
    return x_rot, y_rot
end

function plot_ellipse(a, b, theta)
    t = range(0, 2π, length=100)
    ellipse_x = a * cos.(t)
    ellipse_y = b * sin.(t)
    ellipse_x_rot = cos(theta) .* ellipse_x - sin(theta) .* ellipse_y
    ellipse_y_rot = sin(theta) .* ellipse_x + cos(theta) .* ellipse_y
    plot!(ellipse_x_rot, ellipse_y_rot, lw=2, lc=:gray, label=false)
    
    semi_major_x = [0, a * cos(theta)]
    semi_major_y = [0, a * sin(theta)]
    plot!(semi_major_x, semi_major_y, lw=2, lc=:orange)

    semi_minor_x = [0, -b * sin(theta)]
    semi_minor_y = [0, b * cos(theta)]
    plot!(semi_minor_x, semi_minor_y, lw=2, lc=:purple)
    
    annotate!([(.4, .38, text("a", :orange, 10, :left)),
               (-.1, .3, text("b", :purple, 10, :left)),
               (0.5 * semi_major_x[2], 0.2 * semi_major_y[2], text("θ", :black, 10, :left))])
end

@gif for t in t_values
    x, y = get_position(t)
    
    p = plot([-2, 2], [0, 0], lw=2, lc=:black, arrow=:arrow, size=(500, 500))
    plot!([0, 0], [-2, 2], lw=2, lc=:black, arrow=:arrow)
    
    plot_ellipse(amplitude_x, amplitude_y, theta)
    
    scatter!([x], [y], markersize=10, legend=false)
    
    plot!([x, 0], [y, y], lw=2, lc=:red)
    plot!([x, x], [y, 0], lw=2, lc=:blue)
    
    scatter!([x], [0], markersize=1, color=:red, legend=false)
    scatter!([0], [y], markersize=1, color=:blue, legend=false)
    
    xlims!(-1, 1)
    ylims!(-1, 1)
    
    annotate!([
        (x+.05, .05, text("B", :blue, 10, :left)),
        (-.05, y+.05, text("A", :red, 10, :right))
    ])
end
