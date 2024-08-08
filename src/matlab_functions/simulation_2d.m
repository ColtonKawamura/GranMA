function simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    %% Molecular Dynamics Simulator (Adapted from Mark D. Shattuck, CCNY)
    % Example command: simulation_2d(100, 1, 1, 1.28, 5000, 5000, 0.01, 5, 1)
    
    % Set up initial conditions and visualization
    % Add random initial velocities
    % Replace periodic boundaries with fixed walls
    % Replace Euler with Velocity Verlet
    % Add "Nplotskip" so every frame is not plotted
    % Add gravity
    %%%%% CHANGE FORCE-DETECTION TO SPEED UP %%%%%
    % Add dissipation during collisions
    % Add dt calculation based on sqrt(m/k)
    % Add Ek(nt) storage inside loop
    
    
    % % % Script Variables for debugging
    % K = 100;
    % M = 1;
    % Bv = .1;
    % w_D = 0.36; % 
    % N = 5000;
    % P = 0.001; % 0.021544 0.046416
    % W = 5;
    % seed = 1;
    
    
    % Create the packing name with the exact number format for P
    packing_name = sprintf('2D_N%d_P%s_Width%d_Seed%d', N, num2str(P), W, seed);
    
    % Create the filename
    % filename_output = sprintf('%s_K%d_Bv%d_wD%d_M%d.mat', packing_name, K, Bv, w_D, M);
    filename_output = sprintf('2D_%s_K%d_Bv%d_wD%.2f_M%d.mat', packing_name, K, Bv, w_D, M);
    
    % if exist(filename_output)
    %     return
    % end
    input_pressure = P;
    load(['in/' packing_name '.mat']);
    
    % Nt = round(.9*Lx/(pi*sqrt(M/K)*0.05))
    
    % N = length(Dn);
    flag = true;
    Lx0 = Lx;
    B=0;
    A = P_target/100;
    
    
    dt = pi*sqrt(M/K)*0.05; %  was pi*sqrt(M/K)*0.05
    c_0 = min(Dn).*sqrt(K/M);
    Nt = round(.9.*(Lx ./ c_0)./(dt))
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
    
    
    % identify wall particles
    left_wall_list = (x<Dn/2);
    right_wall_list = (x>Lx-Dn/2);
    bulk_list = ~(left_wall_list | right_wall_list);
    
    %% Main Loop
    % P = 0; % CAK not sure why this was set to zero....
    for nt = 1:Nt
    
    
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
    [fitted_attenuation, wavenumber, attenuation_fit_line, initial_distance_from_oscillation_output, amplitude_vector, unwrapped_phase_vector, cleaned_particle_index] = ...
        process_gm_fft(driving_amplitude, time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, initial_distance_from_oscillation)
    
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
    [fitted_attenuation, wavenumber, attenuation_fit_line, initial_distance_from_oscillation_output, amplitude_vector, unwrapped_phase_vector, cleaned_particle_index] = ...
        process_gm_fft(driving_amplitude, time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, initial_distance_from_oscillation)
    
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
    
    % plot_white_paper(initial_distance_from_oscillation_output_y_fft, initial_distance_from_oscillation_output_x_fft, amplitude_vector_x, amplitude_vector_y, unwrapped_phase_vector_x, unwrapped_phase_vector_y, cleaned_particle_index_x, cleaned_particle_index_y, x_all, y_all)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Figure of one particle's motion, just for poster purposes
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting
    figure
    hold on
    
    % Plot centered on zero
    position_particles = x_all;
    index_particle_to_plot = 32;
    plot(time_vector, position_particles(index_particle_to_plot,:) - mean(position_particles(index_particle_to_plot,:)))
    position_particles = y_all;
    plot(time_vector, position_particles(index_particle_to_plot,:) - mean(position_particles(index_particle_to_plot,:)))
    
    % Title and axis labels in LaTeX format
    title('Particle Position Over Time', 'Interpreter', 'latex')
    xlabel('Time (s)', 'Interpreter', 'latex')
    ylabel('Position', 'Interpreter', 'latex')
    legend('$\hat{k}$', '$\hat{k}_\perp$' ,'Interpreter', 'latex')
    
    % Adding legend
    legend show
    
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ellipse Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tvec = (1:Nt)*dt;
    [~,isort] = sort(x0);
    ilist = [1,50,105,248,500,1000];
    % figure(1044) 
    tiledlayout(2,3,'TileSpacing','tight')
    iskip = 10;
    
    % What is intent of this loop? Just debugging for particles in ilist?
    for nn = isort(ilist)
        x_temp = x_all(nn,:);
        y_temp = y_all(nn,:);
    
        % Find the first index i0 where x_temp exceeds the midpoint between x0(nn) and its maximum value.
        i0 = find(x_temp>x0(nn)+0.5*(max(x_temp)-x0(nn)),1,'first');
    
        % Calculate the amplitude as the distance from the initial position.
        amp = sqrt((x_temp-x0(nn)).^2+(y_temp-y0(nn)).^2);
    
        % Calculate the angle from the initial position.
        ang = atan2((y_temp-y0(nn)),(x_temp-x0(nn)));
    
        % Shorten the time vector from the found index i0 to the end.
        t_temp = tvec(i0:end);
    
        % Compute the median amplitude.
        medamp = median(amp);
    
        % Create a time array for plotting.
        t = 1:length(amp);
    end
    
    tvec = (1:Nt)*dt;
    omega_D = w_D;
    [~,isort] = sort(x0);
    iskip = 1;
    list = [];
    b_start = 0;
    
    ellipse_stats = zeros(length(isort(1:iskip:end)), 4);
    nn_list = isort(1:iskip:N);
    
    random_plot_sample = randsample(nn_list, 10);
    debug_counter = 1;
    for nn_counter = 1:length(nn_list)
        % waitbar(nn_counter/length(nn_list)); % COMMENTED OUT BECAUSE: On hamming waitbar not supported under the -nojvm startup option.
        nn = nn_list(nn_counter);
        if(~left_wall_list(nn))
            x_temp = x_all(nn,:);
            y_temp = y_all(nn,:);
            if length(unique(x_temp))>100
    
                % find index when particles rises to .5 its max amplitude
                i0 = find(x_temp>x0(nn)+0.5*(max(x_temp)-x0(nn)),1,'first');
    
                % Calculate the amplitude as the Euclidean distance from the initial position (x0(nn), y0(nn)).
                amp = sqrt((x_temp-x0(nn)).^2+(y_temp-y0(nn)).^2);
    
                % Calculate the angle from the initial position to each point in Cartesian coordinates.
                ang = atan2((y_temp-y0(nn)),(x_temp-x0(nn)));
                t_temp = tvec(i0:end);
    
                medamp = median(amp);
                t = 1:length(amp);
    
                % Normalize the x and y coordinates by subtracting the initial position and scaling by A.
                X = (x_temp(100:end)-x0(nn))/A; % why does this start at 100? Should it start at t_temp then plot t_temp, X? Aka time we mark as the rise above .5 max amplitude? Otherise you get a long period of no oscillation
                Y = (y_temp(100:end)-y0(nn))/A;
    
                % Why are these set to their same values? 
                X = X; 
                Y = Y;
    
                % Compute the radius squared for each point relative to the origin.
                R = X.^2+Y.^2;
    
                % Find the maximum radius squared from the array R. This should become the semi-major axis
                R_max = max(R);
    
                % Find the index of the first occurrence where R equals R_max.
                max_ind = find(R == R_max,1,"first");
    
                % Extract the X and Y coordinates corresponding to the maximum radius.
                X_max = X(max_ind);
                Y_max = Y(max_ind);
                af = [];
    
                % *** NEEDS fit_ellipse FUNCTION *** Fit an ellipse to the latter half of the data points, returning a structure with ellipse parameters.
                ellipse_t = fit_ellipse(X(round(end/2):end)',Y(round(end/2):end)');
                % test_x = (x_temp(i0:end)-x0(nn))/A;
                % test_y = (y_temp(i0:end)-y0(nn))/A;
                % ellipse_t = fit_ellipse(test_x',test_y');
                % nn_counter
    
                % Check if the ellipse fitting structure 'ellipse_t' is not empty and contains  parameter 'a' that comes from the fit_ellipse function
                if ~isempty(ellipse_t) && ~isempty(ellipse_t.a)
                    % if nn_counter==336
                    
                    % end
    
                    % Assign the semi-major axis, semi-minor axis, and negative rotation angle from 'ellipse_t' to 'af'
                    af(1) = ellipse_t.a; % sub axis (radius) of the X axis of the non-tilt ellipse
                    af(2) = ellipse_t.b; % sub axis (radius) of the Y axis of the non-tilt ellipse
                    af(3) = -ellipse_t.phi; % orientation in radians of the ellipse (tilt)
                else
                    % Define initial parameter estimates for fitting based on maximum radius and the angle at maximum point
                    a0 = [R_max R_max atan(Y_max/X_max)];
                    c = [0 0]; % Center of transformation set to origin
    
                    % Define the function 'f' to minimize an ellipse fit error function parameterized by 'a'
                    f = @(a) ((((X-c(1))*cos(a(3))+(Y-c(2))*sin(a(3))).^2)/a(1).^2 + (((X-c(1))*sin(a(3))-(Y-c(2))*cos(a(3))).^2)/a(2).^2 -1);
                    %[af fval] = fminsearch(f,a0);
    
                    % Use nonlinear least squares to optimize parameters 'a' starting from 'a0'
                    af = lsqnonlin(f, a0);%, [], [], options);
                end
    
                if ~isempty(af)

                    % Make sure to consider only the abs of the ellipse. for the next step
                    af(1:2) = abs(af(1:2));
    
                    % If the semi-major axis (af(1)) is less than the semi-minor axis (af(2))
                    if af(1)<af(2)
    
                        % Rotate the ellipse by -pi/2 to correct the orientation
                        af(3) = af(3)-pi/2;
    
                        % Swap the semi-major and semi-minor axes to maintain convention
                        temp = af(2);
                        af(2) = af(1);
                        af(1) = temp;
                    end
    
                    % Normalize the rotation angle af(3) within the range -pi/2 to pi/2
                    af(3) = (mod(af(3)+pi/2,pi)-pi/2);
    
    
                    % Store the fitted ellipse parameters (semi-major axis, semi-minor axis, rotation angle, initial position from oscillation) in the ellipse_stats matrix for the current particle.
                    ellipse_stats(nn_counter,:) = [af, x0(nn)];

                    debug_aspect_ratio = af(2)./af(1)

                    % if ismember(nn_counter, random_plot_sample)
                    if debug_aspect_ratio > 1
                        debug_counter = debug_counter +1
                        % if debug_counter>5
                        %     return
                        % end
                        figure(nn_counter), clf, plot((x_temp(100:20:end)-x0(nn))/A,(y_temp(100:20:end)-y0(nn))/A,'ro'), axis equal, hold on
                        % plot(X_max,Y_max,'ro','markerfacecolor','r','markersize',20)
                        debug_counter = debug_counter +1
                        % plot(c(1), c(2), 'r*')
                        axis equal
                        theta_rot = af(3);
                        v = linspace(0, 1);
                        rx = af(1);
                        ry = af(2);
                        x = rx*cos(2*pi*v);
                        y = ry*sin(2*pi*v);
                        xx = x*cos(theta_rot) - y*sin(theta_rot);
                        yy = x*sin(theta_rot) + y*cos(theta_rot);
                        % figure
                        plot(xx, yy,'k-','linewidth',2)
                        plot([0 af(1)*cos(af(3))],[0 af(1)*sin(af(3))],'k-')
                        plot([0 -af(2)*sin(af(3))],[0 af(2)*cos(af(3))],'k-')
                        plot([0 af(1)*cos(af(3))],[0 0],'k--')
                        text(0.4, 0.1, '$a$','Interpreter','latex','fontsize',14)
                        text(-0.1, 0.02, '$b$','Interpreter','latex','fontsize',14)
                        text(0.8,0, '$\theta$','Interpreter','latex','fontsize',14)
                        
                        xlim([1.05*min((x_temp(100:20:end)-x0(nn))/A), 1.05*max((x_temp(100:20:end)-x0(nn))/A)])
                        ylim([1.05*min((y_temp(100:20:end)-y0(nn))/A), 1.05*max((y_temp(100:20:end)-y0(nn))/A)])
                        xlabel('$\Delta x / A$','Interpreter','latex','fontsize',18)
                        ylabel('$\Delta y / A$','Interpreter','latex','fontsize',18)
                        title(['$\theta=$' num2str(af(3))], "Interpreter", "latex")
                        % grid
                        % axis([-1  1    -1  1])
                        % axis('equal')
                        % t=0:pi/10:2*pi;
                        % plot(c(1) + af(1)*cos(t), c(2) + af(2)*sin(t), 'r')
                        % axis([-2 2 -2 2])

                    end
                end
            end
        end
    end
    
    % Copy the ellipse_stats matrix to a new variable for processing.
    ellipse_stats_nonzero = ellipse_stats; % (semi-major axis, semi-minor axis, rotation angle, initial position from oscillation)
    
    % Convert the rotation angles from radians to degrees for easier interpretation.
    ellipse_stats_nonzero(:,3) = ellipse_stats_nonzero(:,3)*180/pi;
    
    % Ensure the semi-major and semi-minor axes values are positive.
    ellipse_stats_nonzero(:,1:2) = abs(ellipse_stats_nonzero(:,1:2));
    
    % Remove any rows where the semi-major axis is zero.
    ellipse_stats_nonzero = ellipse_stats_nonzero(ellipse_stats_nonzero(:,1)~=0,:);
    
    % Calculate the histogram of aspect ratios (semi-minor axis divided by semi-major axis) using bins of 0.05.
    [asp_rat_counts,asp_rat_bins] = histcounts((ellipse_stats_nonzero(:,2))./(ellipse_stats_nonzero(:,1)),0:.05:1);
    
    % Calculate the histogram of absolute rotation angles using bins of 5 degrees up to 90 degrees.
    [rot_ang_counts,rot_ang_bins] = histcounts(abs(ellipse_stats_nonzero(:,3)),0:5:90);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Procesces Ellipse Statistics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the limit of what's cleaned data
    
    max_distance_from_wall = round(max(initial_distance_from_oscillation_output_x_fft));
    
    [binned_semi_minor_values, binned_semi_major_values, binned_rot_angle_values, bin_centers] = process_ellipse_data(ellipse_stats_nonzero);
    
    [mean_aspect_ratio, mean_rotation_angles] =plot_ellipse_statistics(binned_semi_minor_values, binned_semi_major_values, binned_rot_angle_values, bin_centers, max_distance_from_wall)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Let's make everything Dimensionless
    diameter_average = 1; % This is only here until all the old packings are updated to have this as an output CAO 1JUN2024
    mass_particle_average = 1;
    attenuation_x_dimensionless = attenuation_x*diameter_average;
    attenuation_y_dimensionless = attenuation_y*diameter_average;
    wavenumber_x_dimensionless = wavenumber_x*diameter_average;
    wavenumber_y_dimensionless = wavenumber_y*diameter_average;
    driving_angular_frequency_dimensionless = w_D*sqrt(mass_particle_average/K);
    gamma_dimensionless = Bv/sqrt(K*mass_particle_average);
    pressure_dimensionless = P;
    % Save the file
    save(['out/simulation_2d/bi_width_effect/' filename_output], 'gamma_dimensionless','ellipse_stats_nonzero', 'asp_rat_bins', 'asp_rat_counts', ...
        'rot_ang_bins', 'rot_ang_counts', 'time_vector', 'index_particles', 'attenuation_x_dimensionless', ...
        'attenuation_y_dimensionless', 'wavenumber_x_dimensionless', 'wavenumber_y_dimensionless', 'wavespeed_x', ...
         'wavespeed_y', 'driving_angular_frequency_dimensionless', 'attenuation_fit_line_x', ...
            'initial_distance_from_oscillation_output_x_fft', 'initial_distance_from_oscillation_output_y_fft', ...
             'amplitude_vector_x', 'amplitude_vector_y', "pressure_dimensionless", "seed", "mean_aspect_ratio", "mean_rotation_angles", "seed", "input_pressure")