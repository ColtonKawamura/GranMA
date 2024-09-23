% function simulation2dGausian(K, M, Bv, w_D, N, P, W, seed)
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
    addpath('./src/matlab_functions') 
    
    % % Script Variables for debugging
    K = 100;
    M = 1;
    Bv = 0;
    w_D = .3 % Low wend is .2 (before hitting wall @ Nt = 20K) high is 1 @ 5000, 2tracking @ omega = .8, P.1
    N = 5000;
    P = 0.001; % 0.021544 0.046416
    W = 5;
    seed = 1;

    save_interval = 10;
    xOut = [];
    yOut = [];
    timeVector = [];


    % Set up Gaussian Envelope for Pulse
    tau = 1/(w_D/(2*pi) )* 2; % ten cycles long, was 5
    sigma = tau  / sqrt(2*log(2)); % spread of the pulse
    t_max = 4*tau; % time when the center of the pulse occurs

    f = @(t,tmax,sigma) exp(-(t-tmax).^2./(sigma.^2)); 

    
    % Create the packing name with the exact number format for P
    packing_name = sprintf('2D_N%d_P%s_Width%d_Seed%d', N, num2str(P), W, seed);
    
    % Create the filename
    % filename_output = sprintf('%s_K%d_Bv%d_wD%d_M%d.mat', packing_name, K, Bv, w_D, M);
    filename_output = sprintf('%s_K%d_Bv%d_wD%.2f_M%d.mat', packing_name, K, Bv, w_D, M);
    
    % if exist(char("out/simulation_2d/K100_no_last_quarter/" + filename_output), 'file')
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
    % Nt = round(.9.*(Lx ./ c_0)./(dt));
    Nt = 35000;
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
    
    %shenanigans
    [x0_sorted,idx]=sort(x0);
    outAmp = [];
    outXinit = [];

    for nt = 1:Nt
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Debug Plotting %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % if mod(nt,100)==0
        %     smoothAmp = smooth(x(idx) - x0(idx), 150, 'sgolay');
        %     minPeakWidth = 1/w_D * .6;
        %     [pks,locs]=findpeaks(smoothAmp, "MinPeakWidth", minPeakWidth);
        %     valid_peaks_idx = find(pks > A * 0.2);
        %     [~, first_peak_idx] = max(locs(valid_peaks_idx));  % find the first peak thats greater than 1/2 amp
        %     firstPeak = pks(valid_peaks_idx(first_peak_idx)); % value of the first peak, store in vector 
        %     peak_index = locs(valid_peaks_idx(first_peak_idx));
        %     if ~isempty(valid_peaks_idx)
        %         % plot(x0(idx), x(idx) - x0(idx), '.', x0(idx), smoothAmp, 'r-', x0(idx(peak_index)), x(idx(peak_index)) - x0(idx(peak_index)), 'o', 'MarkerFaceColor', 'r')
        %         plot(x0(idx), x(idx) - x0(idx), '.', x0(idx), smoothAmp, 'r-', x0(idx(peak_index)), firstPeak, 'o', 'MarkerFaceColor', 'r')
        %         ylim(1.2*[-A,A])
        %         outAmp = [outAmp, firstPeak];
        %         outXinit = [outXinit, x0(idx(peak_index))];
        %     end
        %     drawnow
        % end
        % if mod(nt,10)==0
        %     smoothAmp = smooth(x(idx) - x0(idx), 150, 'sgolay');
        %     minPeakWidth = 1/w_D * .6;
        %     [pks,locs]=findpeaks(smoothAmp, "MinPeakWidth", minPeakWidth);
        %     valid_peaks_idx = find(pks > A * 0.1);
        %     if ~isempty(valid_peaks_idx) && length(valid_peaks_idx) >2
        %         second_peak_idx = valid_peaks_idx(end-1);
        %         second_peak_Amp = pks(second_peak_idx);
        %         second_peak_xInit = x0(idx(locs(second_peak_idx)));
        %         % plot(x0(idx), x(idx) - x0(idx), '.', x0(idx), smoothAmp, 'r-', x0(idx(peak_index)), x(idx(peak_index)) - x0(idx(peak_index)), 'o', 'MarkerFaceColor', 'r')
        %         plot(x0(idx), x(idx) - x0(idx), '.', x0(idx), smoothAmp, 'r-', second_peak_xInit, second_peak_Amp, 'o', 'MarkerFaceColor', 'r')
        %         ylim(1.2*[-A,A])
        %         outAmp = [outAmp, second_peak_Amp];
        %         outXinit = [outXinit, second_peak_xInit];
        %         if second_peak_Amp < outAmp(1) * .7
        %             display("peak dropped")
        %             break
        %         end
        %     end
        %     drawnow
        % end
        [breakOut, outAmp, outXinit]= getAmps(nt, x, x0, idx, w_D, A, outAmp, outXinit);
        % outAmp = [outAmp, second_peak_Amp];
        % outXinit = [outXinit, second_peak_xInit];
        if breakOut
            break; % Break the loop if shouldBreak is true
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% First step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        x_all(:,nt) = x;
        y_all(:,nt) = y;
    
        x  =  x+vx*dt+ax_old.*dt.^2/2;
        y  =  y+vy*dt+ay_old.*dt.^2/2;
    
        x(left_wall_list) = x0(left_wall_list)+A*cos(w_D*((nt)*dt-t_max))*f(nt*dt,t_max,sigma);
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
        
        if mod(nt, save_interval) == 0
            xOut = [xOut, x'];
            yOut = [yOut, y'];
            timeVector = [timeVector, nt*dt];
        end

    end
    
    %%%%% Post Processing %%%%%
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % X Direction Post Processing
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Add the path to the "functions" directory
    addpath('./src/matlab_functions')

    % Dn = Dn';

    % Sorting based on the initial x positions (first column of xOut)
    [~, sortIdx] = sort(xOut(:, 1));  % Get the sorting indices based on initial x positions

    % xOut = xOut(sortIdx, :);
    % yOut = yOut(sortIdx, :);
    meanDiameter = mean(Dn);
    figure
    semilogy(outXinit, outAmp, 'o')

     % Save the file
     save(['out/simulation_2d/' filename_output],  'xOut', 'yOut', 'Dn', 'timeVector', 'tau', 'w_D', 'Bv', 'K', 'M', 'P', 'W', 'seed', 'P_target')


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % Figure of one particle's motion, just for poster purposes
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Vector of target initial distances
    % % target_distances = [11.2249, 62.3078, 128.511, 199.605]; %  simulation_2d(100, 1, .0000000001, 1, 5000, .001, 5, 1)
    % % target_distances = [21.5772, 221.969, 399.241, 599.475]; %  simulation_2d(100, 1, .0000000001, 1, 5000, .1, 5, 1)
    % target_distances = [25, 125, 225, 325]; %  simulation_2d(100, 1, .0000000001, 1, 5000, .1, 5, 1)

    % % Find the indices of the particles closest to the target distances
    % indices_to_plot = zeros(1, length(target_distances)); % Preallocate an array to store the indices
    % for i = 1:length(target_distances)
    %     [~, indices_to_plot(i)] = min(abs(x0 - target_distances(i)));
    % end

    % % Calculate the maximum y-axis value for consistent scaling
    % max_y_value = 0; % Initialize
    % for i = 1:length(indices_to_plot)
    %     index_particle_to_plot = indices_to_plot(i);
        
    %     % Get max values for x and y positions
    %     max_y_value = max(max_y_value, max(abs(x_all(index_particle_to_plot, :) - mean(x_all(index_particle_to_plot, :)))));
    %     max_y_value = max(max_y_value, max(abs(y_all(index_particle_to_plot, :) - mean(y_all(index_particle_to_plot, :)))));
    % end

    % % Set up tiled layout for the plots
    % tiledlayout(length(target_distances),1); % Create tiled layout with rows equal to the number of target distances

    % % Loop over each selected particle to plot its position over time in separate tiles
    % for i = 1:length(indices_to_plot)
    %     index_particle_to_plot = indices_to_plot(i);
        
    %     % Create a tile for each particle
    %     nexttile;
        
    %     % Plot x position
    %     position_particles_x = x_all;
    %     plot(time_vector, position_particles_x(index_particle_to_plot, :) - mean(position_particles_x(index_particle_to_plot, :)), ...
    %         'DisplayName', sprintf('Distance = %.3f (x)', x0(index_particle_to_plot)));
    %     hold on;
        
    %     % Plot y position
    %     position_particles_y = y_all;
    %     plot(time_vector, position_particles_y(index_particle_to_plot, :) - mean(position_particles_y(index_particle_to_plot, :)), ...
    %         'DisplayName', sprintf('(y)'));
        
    %     % Set consistent y-axis limits
    %     ylim([-max_y_value, max_y_value]);
        
    %     % Box, labels, and legends for each tile
    %     box on;
    %     xlabel('Time (s)', 'Interpreter', 'latex');
    %     ylabel('Position', 'Interpreter', 'latex');
    %     legend show;
        
    %     hold off;
    % end
    