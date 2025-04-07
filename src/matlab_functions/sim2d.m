function sim2d(K, M, Bv, w_D, N, P, W, seed, in_path, out_path)
    % Example command: sim2d(100, 1, 1, 1.28, 80000, 0.1, 40, 1, "in/2d_poly_20by20/", "out/junk_yard")
    % sim2d(100, 1, 1, 1, 80000, 0.1, 40, 1, "in/2d_tiled_2000by40/", "out/junk_yard")

    % Create the packing name with the exact number format for P
    packing_name = string(sprintf("2D_N%d_P%s_Width%d_Seed%d", N, num2str(P), W, seed));
    
    % Create the filename
    % filename_output = sprintf('%s_K%d_Bv%d_wD%d_M%d.mat', packing_name, K, Bv, w_D, M);
    filename_output = string(sprintf("%s_K%d_Bv%d_wD%.2f_M%d.mat", packing_name, K, Bv, w_D, M));
    
    % if exist(char("out/simulation_2d/K100_everything_smaller_dt/" + filename_output), 'file')
    %     return
    % end
    input_pressure = P;
    filename = in_path + packing_name + ".mat";  % Concatenate path and filename

    load(filename);
    % load(['in/' packing_name '.mat']);
    
    % Nt = round(.9*Lx/(pi*sqrt(M/K)*0.05))
    
    % N = length(Dn);
    flag = true;
    Lx0 = Lx;
    B=0;
    A = P_target/100;
    
    
    dt = pi*sqrt(M/K)*0.05; %  was pi*sqrt(M/K)*0.05
    c_0 = min(Dn).*sqrt(K/M);
    Nt = round(.9.*(Lx ./ c_0)./(dt));
    ax_old = 0*x;
    ay_old = 0*y;
    vx = 0*x;
    vy = 0*y;
    Ek = zeros(1,Nt);
    Ep = zeros(1,Nt);
    g = 0;
    
    %% initial positions
    
    x0 = x;
    y0 = y;
    
    %% Make neighbor lists with initial spring lengths
    
    skin = 0;
    Zn_list = [];
    neighbor_list_all = [];
    spring_list_all = [];
    for nn = 1:N
        neighbor_list_nn = [];
        spring_list_nn = [];
        for mm = [1:nn-1,nn+1:N]
            dy = y(mm)-y(nn);
            dy = dy - round(dy/Ly)*Ly;
            Dnm = (1+skin)*(Dn(nn) + Dn(mm))/2;
            if(abs(dy) <= Dnm)
                dx = x(mm)-x(nn);
                dnm = dx.^2+dy.^2;
                if(dnm < Dnm^2)
    
                    neighbor_list_nn = [neighbor_list_nn, mm];
                    spring_list_nn = [spring_list_nn, sqrt(dnm)];
    
                end
            end
        end
        neighbor_list_all{nn} = neighbor_list_nn;
        spring_list_all{nn} = spring_list_nn;
        Zn_list = [Zn_list;length(spring_list_nn)];
    end
    
    % [Hessian, eigen_values, eigen_vectors] = HessYale(x, y, Dn, N, Ly, K);
    % sort eigen_vector's columns in ascending order, then short eignvectors accordingly
    % [ascending_eigen_values, eigen_idx] = sort(eigen_values, 1);
    % ascending_eigen_vectors = zeros(size(ascending_eigen_values));
    % for col = 1:size(ascending_eigen_vectors, 2) % Go throuch each column 
    %     ascending_eigen_vectors(:, col) = eigen_vectors(eigen_idx(:, col), col); % grab all the rows from each column in eigen_vectors and 
    % end

    % mode_to_plot = 2; % should just be the column
    % plotEigenmode(x0', y0', eigen_vectors, mode_to_plot)
    % identify wall particles
    left_wall_list = (x<Dn/2);
    right_wall_list = (x>Lx-Dn/2);
    bulk_list = ~(left_wall_list | right_wall_list);
    
    %% Main Loop
    % P = 0; % CAK not sure why this was set to zero....
    [~, idx] = sort(x0);
    for nt = 1:Nt
        % visualizeSim(10000, x, x0, y, y0, idx, A)
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% First step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        x_all(:,nt) = x;
        y_all(:,nt) = y;
    
        x  =  x+vx*dt+ax_old.*dt.^2/2;
        y  =  y+vy*dt+ay_old.*dt.^2/2;
    
        x(left_wall_list) = x0(left_wall_list)+A*sin(w_D*dt*nt);
        y(left_wall_list) = y0(left_wall_list);
        x(right_wall_list) = x0(right_wall_list);
        y(right_wall_list) = y0(right_wall_list);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Interaction detector and Force Law %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        Fx = zeros(1,N);
        Fy = zeros(1,N);
        Zn = zeros(1,N);
    
        for nn = 1:N
            spring_list = spring_list_all{nn};
            neighbor_list = neighbor_list_all{nn};
            for mm_counter = 1:length(neighbor_list)
                mm = neighbor_list(mm_counter);
                dy = y(mm)-y(nn);
                dy = dy - round(dy/Ly)*Ly;
                Dnm = spring_list(mm_counter);
                %             if(abs(dy) <= Dnm)
                dx = x(mm)-x(nn);
                dnm = dx.^2+dy.^2;
                %                 if(dnm < Dnm^2)
                dnm = sqrt(dnm);
    
                F = -K*(Dnm/dnm-1);
    
                dvx = vx(nn)-vx(mm);
                dvy = vy(nn)-vy(mm);
    
                Fx(nn) = Fx(nn)+F.*dx-Bv*dvx;  % particle-particle Force Law
                %                     Fx(mm) = Fx(mm)-F.*dx+Fdiss.*dx/dnm;
                Fy(nn) = Fy(nn)+F.*dy-Bv*dvy;
                %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
                Zn(nn) = Zn(nn) + 1;
                %                     Zn(mm) = Zn(mm) + 1;
                Ep(nt) = Ep(nt) + 0.5*K*(Dnm-dnm)^2;
                %                 end
                %             end
            end
        end
        %
        Fx = Fx - B.*vx;
        Fy = Fy - B.*vy;
    
        Fx(left_wall_list) = 0;
        Fx(right_wall_list) = 0;
    
        Ek(nt) = 1/2*M*sum((vx).^2+(vy).^2);
        Ek(nt) = Ek(nt)/N;
        Ep(nt) = Ep(nt)/N;
    
        ax = Fx./M;
        ay = Fy./M-g;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Second step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        vx = vx+(ax_old+ax).*dt/2;
        vy = vy+(ay_old+ay).*dt/2;
    
        ax_old = ax;
        ay_old = ay;
    end
    
    %%%%% Post Processing %%%%%
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % X Direction Post Processing
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Add the path to the "functions" directory
    addpath('./src/matlab_functions')
    
    % Convert simulation variables to meet function convention
    time_vector = (1:Nt)*dt;
    [~,index_particles] = sort(x0);
    index_oscillating_wall = left_wall_list;
    driving_frequency = w_D/6.2832;
    driving_amplitude=A;
    position_particles = x_all;
    initial_distance_from_oscillation = x0;
    
    % Perform fft fitting
    [fitted_attenuation, wavenumber, attenuation_fit_line, initial_distance_from_oscillation_output, amplitude_vector, unwrapped_phase_vector, cleaned_particle_index, initial_position_y_out, initial_position_z_out] = ...
        process_gm_fft(driving_amplitude, time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, initial_distance_from_oscillation, y0, y0);
    
    % Don't waste time if we didn't get enough data
    if isempty(cleaned_particle_index)
    
        fprintf('Simulation P=%d, Omega=%d, Gamma=%d, Seed=%d did not detect attenuation\n', P, w_D, Bv, seed);
    
        % error('Execution stopped: No attenuation detected.');
        return
    end
    
    attenuation_x = fitted_attenuation;
    attenuation_fit_line_x = attenuation_fit_line;
    wavenumber_x = wavenumber;
    unwrapped_phase_vector_x = unwrapped_phase_vector;
    wavespeed_x = driving_frequency*2*pi*sqrt(M/K)/(wavenumber*1);
    initial_distance_from_oscillation_output_x_fft = initial_distance_from_oscillation_output;
    amplitude_vector_x = amplitude_vector;
    cleaned_particle_index_x = cleaned_particle_index;
    
    % process_gm_fft_freq_density(time_vector, index_particles, index_oscillating_wall, driving_amplitude, position_particles, initial_distance_from_oscillation, driving_frequency)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Y Direction Post Processing
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Convert simulation variables to meet function convention
    time_vector = (1:Nt)*dt;
    index_oscillating_wall = left_wall_list;
    driving_amplitude=A;
    [~,index_particles] = sort(y0);
    position_particles = y_all;
    initial_distance_from_oscillation = x0;
    driving_frequency = w_D/6.2832;
    
    % Perform fft fitting
    [fitted_attenuation, wavenumber, attenuation_fit_line, initial_distance_from_oscillation_output, amplitude_vector, unwrapped_phase_vector, cleaned_particle_index, initial_position_y_out, initial_position_z_out] = ...
        process_gm_fft(driving_amplitude, time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, initial_distance_from_oscillation, y0, y0);
    
    attenuation_y = fitted_attenuation;
    attenuation_fit_line_y = attenuation_fit_line;
    wavenumber_y = wavenumber;
    unwrapped_phase_vector_y = unwrapped_phase_vector;
    wavespeed_y = driving_frequency*2*pi*sqrt(M/K)/(wavenumber*1);
    initial_distance_from_oscillation_output_y_fft = initial_distance_from_oscillation_output;
    amplitude_vector_y = amplitude_vector;
    cleaned_particle_index_y = cleaned_particle_index;
    
    %  process_gm_fft_freq_density(time_vector, index_particles, index_oscillating_wall, driving_amplitude, position_particles, initial_distance_from_oscillation)
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Figure of attenuation / wavenumber, just for poster purposes
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wrapped_phase_vector_x = mod(unwrapped_phase_vector_x, 2*pi);
    wrapped_phase_vector_y = mod(unwrapped_phase_vector_y, 2*pi);
    figure;
    scatter(initial_distance_from_oscillation_output_x_fft, wrapped_phase_vector_x, 'o');
    grid on;
    hold on;  % Keep the plot for adding the fitted line
    box on;  
    scatter(initial_distance_from_oscillation_output_y_fft, wrapped_phase_vector_y, 'o');
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Figure of one particle's motion, just for poster purposes
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    target_distance = 21;
    [~, index_to_plot] = min(abs(x0 - target_distance));

    % Plotting
    figure
    
    % Plot centered on zero
    position_particles = x_all;
    index_particle_to_plot = index_to_plot;
    plot(time_vector, position_particles(index_particle_to_plot,:) - mean(position_particles(index_particle_to_plot,:)))
    hold on
    position_particles = y_all;
    plot(time_vector, position_particles(index_particle_to_plot,:) - mean(position_particles(index_particle_to_plot,:)))
    
    % title('Particle Position Over Time', 'Interpreter', 'latex')
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 15)
    ylabel('$ A(x)$', 'Interpreter', 'latex', "Rotation", 0, 'FontSize', 15)
    legend('$A_{||}$', '$A_{\perp}$' ,'Interpreter', 'latex', 'FontSize', 15)
    legend show
    grid on
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % FFT of a single partcile 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    position_particles = x_all; 
    data = position_particles(index_particle_to_plot,:) - mean(position_particles(index_particle_to_plot,:));
    fps = 1 / mean(diff(time_vector));
    marker_color = 'b';
    plotfft(data, fps, marker_color)
    grid on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Figure of one particle's motion, just for poster purposes
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vector of target initial distances
    % target_distances = [11.2249, 62.3078, 128.511, 199.605]; %  simulation_2d(100, 1, .0000000001, 1, 5000, .001, 5, 1)
    target_distances = [21.5772, 221.969, 399.241, 599.475]; %  simulation_2d(100, 1, .0000000001, 1, 5000, .1, 5, 1)

    % Find the indices of the particles closest to the target distances
    indices_to_plot = zeros(1, length(target_distances)); % Preallocate an array to store the indices
    for i = 1:length(target_distances)
        [~, indices_to_plot(i)] = min(abs(x0 - target_distances(i)));
    end

    % Calculate the maximum y-axis value for consistent scaling
    max_y_value = 0; % Initialize
    for i = 1:length(indices_to_plot)
        index_particle_to_plot = indices_to_plot(i);
        
        % Get max values for x and y positions
        max_y_value = max(max_y_value, max(abs(x_all(index_particle_to_plot, :) - mean(x_all(index_particle_to_plot, :)))));
        max_y_value = max(max_y_value, max(abs(y_all(index_particle_to_plot, :) - mean(y_all(index_particle_to_plot, :)))));
    end
    figure
    % Set up tiled layout for the plots
    tiledlayout(length(target_distances),1); % Create tiled layout with rows equal to the number of target distances

    % Loop over each selected particle to plot its position over time in separate tiles
    for i = 1:length(indices_to_plot)
        index_particle_to_plot = indices_to_plot(i);
        
        % Create a tile for each particle
        nexttile;
        
        % Plot x position
        position_particles_x = x_all;
        plot(time_vector, position_particles_x(index_particle_to_plot, :) - mean(position_particles_x(index_particle_to_plot, :)), ...
            'DisplayName', sprintf('Distance = %.3f (x)', x0(index_particle_to_plot)));
        hold on;
        
        % Plot y position
        position_particles_y = y_all;
        plot(time_vector, position_particles_y(index_particle_to_plot, :) - mean(position_particles_y(index_particle_to_plot, :)), ...
            'DisplayName', sprintf('(y)'));
        
        % Set consistent y-axis limits
        ylim([-max_y_value, max_y_value]);
        
        % Box, labels, and legends for each tile
        box on;
        xlabel('Time (s)', 'Interpreter', 'latex');
        ylabel('Position', 'Interpreter', 'latex');
        legend show;
        
        hold off;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Let's make everything Dimensionless
    % diameter_average = 1; % This is only here until all the old packings are updated to have this as an output CAO 1JUN2024
    mass_particle_average = 1;
    attenuation_x_dimensionless = attenuation_x*diameter_average;
    attenuation_y_dimensionless = attenuation_y*diameter_average;
    wavenumber_x_dimensionless = wavenumber_x*diameter_average;
    wavenumber_y_dimensionless = wavenumber_y*diameter_average;
    driving_angular_frequency_dimensionless = w_D*sqrt(mass_particle_average/K);
    gamma_dimensionless = Bv/sqrt(K*mass_particle_average);
    pressure_dimensionless = P;
    % Save the file
    % save(['out/simulation_2d/K100_everything_smaller_dt/' filename_output], 'gamma_dimensionless', 'time_vector', 'index_particles', 'attenuation_x_dimensionless','attenuation_y_dimensionless', 'wavenumber_x_dimensionless', 'wavenumber_y_dimensionless', 'wavespeed_x','wavespeed_y', 'driving_angular_frequency_dimensionless', 'attenuation_fit_line_x', 'initial_distance_from_oscillation_output_x_fft', 'initial_distance_from_oscillation_output_y_fft', 'amplitude_vector_x', 'amplitude_vector_y', "pressure_dimensionless", "seed", "input_pressure", "unwrapped_phase_vector_x", "unwrapped_phase_vector_y")
    save_path = fullfile(out_path, filename_output);

    % save(save_path, 'gamma_dimensionless', 'time_vector', 'index_particles', 'attenuation_x_dimensionless','attenuation_y_dimensionless', 'wavenumber_x_dimensionless', 'wavenumber_y_dimensionless', 'wavespeed_x','wavespeed_y', 'driving_angular_frequency_dimensionless', 'attenuation_fit_line_x', 'initial_distance_from_oscillation_output_x_fft', 'initial_distance_from_oscillation_output_y_fft', 'amplitude_vector_x', 'amplitude_vector_y', "pressure_dimensionless", "seed", "input_pressure", "unwrapped_phase_vector_x", "unwrapped_phase_vector_y")
    
    % Not doing ellispse stats, so save as zeros
    ellipse_stats_nonzero = zeros(1, 6);
    asp_rat_bins = zeros(1, 6);
    asp_rat_counts = zeros(1, 6);
    rot_ang_bins = zeros(1, 6);
    rot_ang_counts = zeros(1, 6);
    mean_aspect_ratio = 0;
    mean_rotation_angles = 0;
    y0_yfft = 0;
    y0_xfft = 0;
    
    save( save_path,'gamma_dimensionless','ellipse_stats_nonzero', 'asp_rat_bins', 'asp_rat_counts', ...
        'rot_ang_bins', 'rot_ang_counts', 'time_vector', 'index_particles', 'attenuation_x_dimensionless', ...
        'attenuation_y_dimensionless', 'wavenumber_x_dimensionless', 'wavenumber_y_dimensionless', 'wavespeed_x', ...
         'wavespeed_y', 'driving_angular_frequency_dimensionless', 'attenuation_fit_line_x', ...
            'initial_distance_from_oscillation_output_x_fft', 'initial_distance_from_oscillation_output_y_fft', ...
             'amplitude_vector_x', 'amplitude_vector_y', "pressure_dimensionless", "seed", "mean_aspect_ratio", "mean_rotation_angles", "seed", "input_pressure", "unwrapped_phase_vector_x", "unwrapped_phase_vector_y", "y0_yfft", "y0_xfft");
    % save(['out/simulation_2d/K100_everything3/' filename_output])