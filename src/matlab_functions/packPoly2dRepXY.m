function packPoly2dRepXY(N, K, D, G, M, P_target, W_factor, seed, plotit, x_mult, y_mult, calc_eig, save_path)
    %Function to create 2D packing with following input parameters:
    % N, Number of Particles that will be repreated 10 times
    % K, spring constant
    % D, Average Diameter
    % G, Ratio of large to small particles (typically 1.4)
    % M, mass of particles
    % P_thres, targeted threshold pressure
    % W_factor, Factor of the width vs number of particles d
    % calc_eig, boolean wheter or not you want to caculate and save eigen stuff 
    % packBi2dRepXY(26,100, 1, 1.5, 1, .1, 5, 1, true, 1, 1, false, 'in/damped_eig_test/')
    
    %% Set up section
    rng(seed)

    N_new = N * x_mult*y_mult; % number of times x repeated times y repeated
    W_new = W_factor * y_mult; % just number of times y is repeated
    % filename = ['in/2D_N' num2str(N_new) '_P' num2str(P_target) '_Width' num2str(W_new) '_Seed' num2str(seed) '.mat']
    filename = [save_path, '2D_N' num2str(N_new) '_P' num2str(P_target) '_Width' num2str(W_new) '_Seed' num2str(seed) '.mat'];

    if exist(filename)
        display("file already exsists!")
        return
    end
    
    Lx = N*D/W_factor; % box width
    Ly = N*D/2; %starting box height
    Bv = 0.1; % dissipation factor
    B = 0.5; % absolute dissipation
    T = 1; % temperature factor
    Nsmall = N/2; % Number of small
    Nbig = N/2; % Number of big
    
    Dn=rand(1,N); %randomize diameter
    [~, i]=sort(Dn);
    for k=1:N
        Dn(i(k))=D-G/2 + G*Dn(i(k)); %starting with the smallest diameter of Small particles, assign a diamter that between D and D+G
    end
    
    
    %% Physical Parameters
    g = 0;
    
    P = 0;
    P_fast_grow = P_target/50;
    
    r = P_target;
    r_fast = 0.01;
    flag = true;
    flag2 = true;
    fast_compress_flag = true;
    
    %% Display Parameters
    % plotit = 1;  % plot ?
    plot_KE = 0;
    Nplotskip = 100;  % number of timesteps to skip before plotting
    Ncellupdate = 1;
    
    %% Simulation Parmeters
    dt = pi*sqrt(M/K)*0.05;
    Nt = 1e8; % Number of steps
    
    %% Initial Conditions
    [x, y] = ndgrid(D/2:D:Lx-D/2,D/2:D:Ly-D/2);
    [~, ii] = sort(rand(1,numel(x)));
    x = x(ii(1:N));
    y = y(ii(1:N));
    
    vx = sqrt(T)*randn(1,N);
    vx = vx - mean(vx);
    vy = sqrt(T)*randn(1,N);
    vy = vy - mean(vy);
    ax_old = 0*x;
    ay_old = 0*y;
    Ek = zeros(1,Nt);
    Ep = zeros(1,Nt);
    
    %% Verlet cell parameters
    raw_cell_width = 2*G*D;
    num_cells = (Lx/raw_cell_width)
    num_cells = round(num_cells)
    cell_width = Lx/num_cells
    jj_max = ceil(Lx/cell_width);
    
    cell_num = ceil(x/cell_width);
    for jj = 1:jj_max
        cell_list{jj} = find(cell_num == jj);
    end
    
    %% Setup Plotting
    if plotit
        figure(1), clf;
        h=zeros(1,2*N);
        for np = 1:N
            h(np) = rectangle('Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
            % h(np+N)=rectangle('Position',[Lx Ly Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
        end
        axis('equal');
        axis([0 Lx 0 Ly]);
        % pause
    
        figure(2), clf;
    %     figure(3), clf;
    end
    %% Main Loop
    for nt = 1:Nt
    
        % plot particles
        if(plotit && mod(nt,Nplotskip) == 0)
            if flag
                figure(1);
                
                % Main loop for updating particle positions
                for np = 1:N
                    set(h(np), 'Position', [x(np) - 0.5 * Dn(np), y(np) - 0.5 * Dn(np), Dn(np), Dn(np)]);  % [x, y, width, height]
                end
                
                Np = N; 
                
                % ii = find(y < Dn/2);  % Bottom boundary condition
                % for nn = 1:length(ii)
                %     np = ii(nn);
                %     set(h(nn + Np), 'Position', [x(np) - 0.5 * Dn(np), y(np) - 0.5 * Dn(np) + Ly, Dn(np), Dn(np)]);
                % end
                % Np = Np + length(ii);
                % 
                % ii = find(y > Ly - Dn/2);  % Top boundary condition
                % for nn = 1:length(ii)
                %     np = ii(nn);
                %     set(h(nn + Np), 'Position', [x(np) - 0.5 * Dn(np), y(np) - 0.5 * Dn(np) - Ly, Dn(np), Dn(np)]);
                % end
                % Np = Np + length(ii);
                % 
                % ii = find(x > Lx - Dn/2);  % Right boundary condition
                % for nn = 1:length(ii)
                %     np = ii(nn);
                %     set(h(nn + Np), 'Position', [x(np) - 0.5 * Dn(np) - Lx, y(np) - 0.5 * Dn(np), Dn(np), Dn(np)]);
                % end
                % 
                % ii = find(x < Dn/2);  % Left boundary condition
                % for nn = 1:length(ii)
                %     np = ii(nn);
                %     set(h(nn + Np), 'Position', [x(np) - 0.5 * Dn(np) + Lx, y(np) - 0.5 * Dn(np), Dn(np), Dn(np)]);
                % end
                
                figure(1);
                ylim([0, Ly]);
                title(num2str(Ly));
                
            end
            
            figure(2), semilogy(nt, Ek(nt-1), 'ro');
            hold on, semilogy(nt, Ep(nt-1), 'bs');
            hold on, plot(nt,P,'kx')
    
    %         figure(3), plot(nt,Ly,'ro',nt,Ly_min,'bx',nt,Ly_max,'kp'), hold on
    
            drawnow;
    
        elseif (plot_KE && mod(nt,Nplotskip) == 0)
            figure(1), plot(x,y,'k.')
            drawnow
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% First step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        x  =  x+vx*dt+ax_old.*dt.^2/2;
        y  =  y+vy*dt+ay_old.*dt.^2/2;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Re-assign particles to cells %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        if flag2 == true || mod(nt,Ncellupdate) == 0
            flag2 = false;
            cell_num = ceil(x/cell_width);
            for jj = 1:jj_max
                cell_list{jj} = find(cell_num == jj);
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Interaction detector and Force Law %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        Fx = zeros(1,N);
        Fy = zeros(1,N);
        Zn = zeros(1,N);
    
        for jj = 1:jj_max
    
            if jj == 1
                mm_list = [cell_list{1} cell_list{2} cell_list{jj_max}];
            elseif jj == jj_max
                mm_list = [cell_list{jj_max-1} cell_list{jj_max} cell_list{1}];
            else
                mm_list = [cell_list{jj-1} cell_list{jj} cell_list{jj+1}];
            end
    
            for nn = cell_list{jj}
    
                mm_list_nn = mm_list(mm_list~=nn);
    
                for mm = mm_list_nn
                    dy = y(mm)-y(nn);
                    dy = dy - round(dy/Ly)*Ly;
                    Dnm = (Dn(nn) + Dn(mm))/2;
                    if(abs(dy) <= Dnm)
                        dx = x(mm)-x(nn);
                        dx = dx - round(dx/Lx)*Lx;
                        dnm = dx.^2+dy.^2;
                        if(dnm < Dnm^2)
                            dnm = sqrt(dnm);
    
                            F = -K*(Dnm/dnm-1);
    
                            m_red = M*M/(M+M);
                            v_dot_r=((vx(nn)-vx(mm))*dx + (vy(nn)-vy(mm))*dy);
                            Fdiss = Bv * m_red * v_dot_r;
    
                            Fx(nn) = Fx(nn)+F.*dx-Fdiss.*dx/dnm;  % particle-particle Force Law
                            Fy(nn) = Fy(nn)+F.*dy-Fdiss.*dy/dnm;
                            Zn(nn) = Zn(nn) + 1;
                            Ep(nt) = Ep(nt) + 0.5*K*(Dnm-dnm)^2;
                        end
                    end
                end
            end
    
        end
    
        Fx = Fx - B.*vx;
        Fy = Fy - B.*vy;
    
        LW_contacts = x<Dn/2;
        % Fx = Fx-K*(x-Dn/2).*(LW_contacts);  % Left wall
        % Fy = Fy-K*(y-D/2).*(y<D/2);  % Bottom wall
    
        RW_contacts = x>Lx-Dn/2;
        % Fx = Fx-K*(x-(Lx-Dn/2)).*(RW_contacts);  % Right wall
        % Fy = Fy-K*(y-(Ly-D/2)).*(y>Ly-D/2);  % Top wall
        
        y=mod(y,Ly); %periodic boundaries for top and bottom
        x=mod(x,Lx); %periodic boundaries for left and right
    
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
        
        no_cont_list = (Zn == 0);% & ~LW_contacts & ~RW_contacts);
        vx(no_cont_list) = 0;
        vy(no_cont_list) = 0;
        ax(no_cont_list) = 0;
        ay(no_cont_list) = 0;
    % 
    
        ax_old = ax;
        ay_old = ay;
    
        tot_contacts = sum(Zn)/2;
        wall_contacts = sum(LW_contacts) + sum(RW_contacts);
        num_rattlers = sum(Zn==0);
        excess_contacts = tot_contacts + wall_contacts - 2*(N-num_rattlers)+1-1;
    
        P = sqrt(2*Ep(nt)/(K));
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPRESSION DECISIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        if fast_compress_flag
            if P<(P_target/50)
                Ly = Ly * (1-r_fast);
                y = y.*(1-r_fast);
                flag = true;
                flag2 = true;
                nt_compress = nt;
            elseif (P<P_target && Ek(nt)<1e-8)
                Ly = Ly * (1-r_fast);
                y = y.*(1-r_fast);
                flag = true;
                flag2 = true;
                nt_compress = nt;
            elseif (P>P_target && Ek(nt)<1e-10)
                Ly = Ly * (1+r_fast);
                y = y.*(1+r_fast);
                flag = true;
                flag2 = true;
                nt_compress = nt;
                fast_compress_flag = false;
            end
        else
            if P<P_fast_grow
                Ly = Ly * (1-r);
                y = y.*(1-r);
                flag = true;
                flag2 = true;
                nt_compress = nt;
            elseif (P<P_target && Ek(nt)<1e-8)
                Ly = Ly * (1-r);
                y = y.*(1-r);
                flag = true;
                flag2 = true;
                nt_compress = nt;
            elseif(P>P_target && Ek(nt)<1e-20)%sum(Cn)/2>(3*sum(Cn>0)-2))
                break;
            end
        end
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    
    N_repeated = x_mult;
    x_repeated = [];
    y_repeated = [];
    Dn_repeated = [];
    
    for i = 0:N_repeated-1
        x_shifted = x + i * Lx;
        y_shifted = y + i * Ly;
        x_repeated = [x_repeated, x_shifted];
        y_repeated = [y_repeated, y];
        Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
    end
    
    N = N * N_repeated;
    x = x_repeated;
    y = y_repeated;
    Dn = Dn_repeated; % Set Dn to the repeated diameters
    
    % Time for y
    N_repeated = y_mult;
    x_repeated = [];
    y_repeated = [];
    Dn_repeated = [];
    
    for i = 0:N_repeated-1
        y_shifted = y + i * Ly;
        x_repeated = [x_repeated, x];
        y_repeated = [y_repeated, y_shifted];
        Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
    end
    
    N = N * N_repeated;
    x = x_repeated;
    y = y_repeated;
    Dn = Dn_repeated; % Set Dn to the repeated diameters
    W_factor = W_factor * N_repeated;
    % figure;
    % hold on;
    % axis equal;
    % for np = 1:N
    %     rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], 'Curvature', [1, 1], 'EdgeColor', 'b');
    % end
    % axis([0, N_repeated * Lx, 0, Ly]);
    % hold off;
    % pause
    figure;
    hold on;
    axis equal;
    axis([0, x_mult * Lx, 0, Ly*y_mult]);
    
    % Loop through each particle and plot its rectangle
    for np = 1:N
        % Compute the position for the rectangle using the center coordinates (x, y)
        % and the diameter Dn (width and height).
        rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], ...
            'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5); % Optional LineWidth for better visibility
    end
    % Update the figure display
    drawnow;
    hold off;
    
    disp(['number of excess contacts = ' num2str(sum(Zn)/2 + sum(LW_contacts) + sum(RW_contacts) - 2*N)])
    Lx = Lx * x_mult; 
    Ly = Ly * y_mult;
    
    if calc_eig == true
        positions = [x',y'];
        radii = Dn./2;
        [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
        Hessian = hess2d(positions, radii, K, Ly, Lx);
        [eigen_vectors, eigen_values ] =  eig(Hessian);
        save(filename, 'x', 'y', 'Dn', 'Lx', 'Ly', 'K', 'P_target', 'P', 'N', 'eigen_vectors', 'eigen_values');
    else
        save(filename, 'x', 'y', 'Dn', 'Lx', 'Ly', 'K', 'P_target', 'P', 'N');
    end
    

