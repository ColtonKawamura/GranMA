using Random
using Plots
using LaTeXStrings
using Statistics
plotlyjs()

function main()
    N = 100
    W = 5
    diameter_average = 1
    diameter_spread = .4
    bi_disperse = true
    K = 100

    particle_diameters, positions = create_packings_3D(N,W, diameter_average, diameter_spread, bi_disperse)
    plot_packing(positions, particle_diameters)
    compression_phase(K, particle_diameters, positions)
end

function compression_phase(K, particle_diameters, positions)

    N = size(positions,1)

    # Initialized the velocities
    temperature = 1
    mass = 1

    velocities = sqrt(temperature/mass).*rand(N,3) # based on equipartition theroem: for one dimension, .5*m*rms(v^2) = .5*k_B*T -> v = sqrt(k_B*T/m)
    velocities .-= mean(velocities, dims=1) ## Want to net temperature of the sytem to be zero, so subtract the mean so they all add to zero simulation; X .-= 2 means X = X-2

    # Initialize the acceleration as zero
    accelerations = velocities.*0

    # Verlet Integration
    N_steps = 1E8
    dt = pi*sqrt(mass/K)*0.05;

    for step_N = 1:N_steps
        positions .+= velocities .* dt .+ 0.5 .* accelerations .* dt^2
    end

end

function create_packings_3D(N,W, diameter_average, diameter_spread, bi_disperse)

    if bi_disperse
        diameter_large = diameter_average + diameter_spread
        diameter_small = diameter_average
    end

    # Calulate the heigh to of the simulation box and scale to the size of the largest particle
    x_limit_upper = W*diameter_large
    y_limit_upper = W*diameter_large
    z_limit_upper = N/(W^2)*diameter_large

    ##  Create a rank 3 tensor that is x_limit_upper by y_limit_upper by z_limit_upper where each point is the x,y,z coordinate that is evenly spaced by diameter_large 
    # Generate ranges for x, y, and z coordinates so they don't extend past the boundary
    x_range = collect(range(diameter_large/2, x_limit_upper-diameter_large/2, step=diameter_large)) 
    y_range = collect(range(diameter_large/2, y_limit_upper-diameter_large/2, step=diameter_large))
    z_range = collect(range(diameter_large/2, z_limit_upper-diameter_large/2, step=diameter_large))

    # Create 3D grid of coordinates
    x_coords = [x for x in x_range, y in y_range, z in z_range] # creates length(z_range) matricies that have dimensions length(x) by length(y) where each element is the x-position
    y_coords = [y for x in x_range, y in y_range, z in z_range] # puts in all the values for x (in x_range) into each element element of y (in y_range) and each one of those into each element of z (in z_range)
    z_coords = [z for x in x_range, y in y_range, z in z_range] 

    # concatanate into a matrix where each row is the x,y,z coordinate of each particle has 
    positions = hcat(reshape(x_coords, :, 1), reshape(y_coords, :, 1), reshape(z_coords, :, 1))

    # Readjust the number of particles to those that actually fit in the box
    N_actual = size(positions,1)
    
    ## Take the new N and then randomly assign a partcile diamter to reach between diameter_large and diameter_small
    particle_diameters = diameter_small .+ (diameter_large - diameter_small) .* rand(N_actual)

    return particle_diameters, positions
end

function plot_packing(positions, particle_diameters)
    # Create the 3D plot
    x = positions[:, 1]
    y = positions[:, 2]
    z = positions[:, 3]
    marker_sizes = particle_diameters * 10  # Adjust the scaling factor if needed for better visualization

    # Readjust the number of particles to those that actually fit in the box
    N_actual = size(positions,1)

    # Generate hover texts including particle number and coordinates
    hover_texts = ["Particle $i: (x, y, z) = ($(x[i]), $(y[i]), $(z[i]))" for i in 1:N_actual]
    plot_before_compression = scatter3d(x, y, z, markersize=marker_sizes, label="", xlabel="X", ylabel="Y", zlabel="Z", title="3D Scatter Plot of Particles", hover=hover_texts)
    display(plot_before_compression)
end
