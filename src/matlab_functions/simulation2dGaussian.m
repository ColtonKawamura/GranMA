function simulation2dGaussian(K, M, Bv, w_D, N, P, W, seed)
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
    %  K = 100;
    %  M = 1;
    %  Bv = 0;
    %  w_D = .2 % Low wend is .2 (before hitting wall @ Nt = 20K) high is 1 @ 5000, 2tracking @ omega = .8, P.1
    %  N = 5000;
    %  P = 0.001; % 0.021544 0.046416
    %  W = 5;
    %  seed = 1;

    % save_interval = 10;
    % xOut = [];
    % yOut = [];
    % timeVector = [];
    nt_out = [];


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
    Nt = 30000;
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
        [breakOut, outAmp, outXinit, nt_out]= getAmps(nt, x, x0, idx, w_D, A, outAmp, outXinit, Lx, nt_out);
        % [breakOut, outAmp, outXinit, nt_out]= getAmpsK(nt, x, x0, idx, w_D, A, outAmp, outXinit, Lx, nt_out); % use this for full motion
        % [breakOut, outAmp, outXinit, nt_out]= getAmpsGIF(nt, x, x0, idx, w_D, A, outAmp, outXinit, Lx, nt_out);
        if any(x(idx(end-200:end)) - x0_sorted(end-200:end) > A * 0.01)
            display("wave reached back wall")
            breakOut = true;
        end

        if breakOut
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% First step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        x_all(:,nt) = x;
        y_all(:,nt) = y;
    
        x  =  x+vx*dt+ax_old.*dt.^2/2;
        y  =  y+vy*dt+ay_old.*dt.^2/2;
    
        x(left_wall_list) = x0(left_wall_list)+A*cos(w_D*((nt)*dt-t_max))*f(nt*dt,t_max,sigma);
        % x(left_wall_list) = x0(left_wall_list)+A*cos(w_D*((nt)*dt-t_max));
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
        
        % This is for exporting time-space vectors
        % if mod(nt, save_interval) == 0
        %     xOut = [xOut, x'];
        %     yOut = [yOut, y'];
        %     timeVector = [timeVector, nt*dt];
        % end

    end
    
    %%%%% Post Processing %%%%%
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % X Direction Post Processing
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Add the path to the "functions" directory
    addpath('./src/matlab_functions')

    % Dn = Dn';

    % Sorting based on the initial x positions (first column of xOut)
    % [~, sortIdx] = sort(xOut(:, 1));  % Get the sorting indices based on initial x positions

    % xOut = xOut(sortIdx, :);
    % yOut = yOut(sortIdx, :);
    meanDiameter = mean(Dn);
    attenuation =  getSlopeLog(outXinit, outAmp);
    % [slope1, slope2] = getSlopeK(outXinit, outAmp);
    % attenuation = min([abs(slope1), abs(slope2)])

    wavespeed = getSlope(nt_out * dt, outXinit) % particle diamters / time
 

    mass = M;
    spring_constant = K;
    omega = w_D * sqrt(mass / spring_constant );
    gamma = Bv / sqrt(spring_constant * mass);
    attenuation = - meanDiameter * attenuation
    width = W;
    pressure = P_target;
    pressure_actual = P;
    wavespeed = wavespeed / (sqrt(spring_constant / mass) * meanDiameter )
    wavenumber = omega / wavespeed;


     save(['out/simulation_2d/gausGetAmps_no_backtrack_moreFreqsV2/' filename_output], 'meanDiameter', 'tau', 'omega', 'gamma', 'spring_constant', 'mass', 'pressure', 'width', 'seed', 'pressure_actual', 'attenuation', 'wavespeed', 'wavenumber', 'dt')

