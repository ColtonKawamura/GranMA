function packBi3dRepXYZ(N, K, D, G, M, P_target, W_factor, seed, plotit, x_mult, y_mult, z_mult, calc_eig, save_path)
% packBi3dRepXYZ - Generates a 3D particle packing arrangement.
%
% Syntax:
%   packBi3dRepXYZ(N, K, D, G, M, P_target, W_factor, seed, plotit, x_mult, y_mult, z_mult, calc_eig, save_path)
%
% Inputs:
%   N        - Number of particles
%   K        - Spring constant
%   D        - Average particle diameter
%   G        - Ratio of large to small particle diameters (typically 1.4)
%   M        - Mass of particles
%   P_target - Targeted threshold pressure
%   W_factor - Factor controlling width relative to number of particles
%              (Assumes a square width for now)
%   seed     - Random seed for reproducibility
%   plotit   - Boolean flag to visualize the packing (true/false)
%   x_mult   - Scaling factor for x-dimension tiling
%   y_mult   - Scaling factor for y-dimension tiling
%   z_mult   - Scaling factor for z-dimension tiling
%   calc_eig - Boolean flag to compute eigenvalues (true/false)
%   save_path- Path to save output data
%
% Example:
%   packBi3dRepXYZ(126, 100, 1, 1.4, 1, 0.1, 5, 1, false, 100, 10, 10, false, 'in/3d/3d_junkyard')
%   % Generates a 5x5x5 tiled packing structure

%% Set up section
rng(seed)

N_new = N * x_mult * y_mult * z_mult;
W_new = W_factor * y_mult; % !! Currently, this only inclues square-width packings (y_mult = z_mult)

% filename = ['in/3D_N' num2str(N) '_P' num2str(P_target) '_Width' num2str(W_factor) '_Seed' num2str(seed) '.mat'];
filename = [save_path, '3D_N' num2str(N_new) '_P' num2str(P_target) '_Width' num2str(W_new) '_Seed' num2str(seed) '.mat'];
filename_tile = [save_path, '3D_tile_N' num2str(N) '_P' num2str(P_target) '_Width' num2str(W_factor) '_Seed' num2str(seed) '.mat'];

if exist(filename)
    display("Output Already in Folder, Stopping.")
    return
end

Lx = N*D/W_factor^2; % box depth - how many layer's of Ly*Lz's can fit
Ly = N*D/2; %starting box height
Lz = W_factor*D;
Bv = 0.1; % dissipation factor
B = 0.5; % absolute dissipation
T = 1; % temperature factor
Nsmall = N/2; % Number of small
Nbig = N/2; % Number of big

Dn=rand(1,N); %randomize diameter
[~, i]=sort(Dn);
for k=1:Nsmall
    Dn(i(k))=D; %starting with the smallest diameter of Small particles, assign a diamter
end
for k=1:Nbig
    Dn(i(k+Nsmall))=D*G; %starting with the smallest diameter of Big particles, assign a diamter
end

%% Physical Parameters
g = 0;

P = 0;
P_fast_grow = P_target/50;

r = P_target;
r_fast = 0.01;
flag = true;
flag2 = true; % flag for re-assinging particles to cells
fast_compress_flag = true; % !! was set to true for debugging

%% Display Parameters
% plotit = 1;  % plot ?
plot_KE = 0;
Nplotskip = 200;  % number of timesteps to skip before plotting
Ncellupdate = 1;
%% Simulation Parmeters
dt = pi*sqrt(M/K)*0.05;
Nt = 1e8; % Number of steps

%% Initial Conditions
[x, y] = ndgrid(D/2:D:Lx-D/2,D/2:D:Ly-D/2); % make matrix of center-positions
[~, ii] = sort(rand(1,numel(x))); % create random uniform numbers for each particle, then get indices !! not sure why not just use randperm(N,N)
x = x(ii(1:N)); % shuffle the x-positions randomly
y = y(ii(1:N)); 
z = Lz*rand(size(x)); % z doesn't need to be a grid - make it uniform random

vx = sqrt(T)*randn(1,N);
vx = vx - mean(vx); % make the velocities random, but shift to conserve momentum (total v =0)
vy = sqrt(T)*randn(1,N);
vy = vy - mean(vy);
vz = sqrt(T)*randn(1,N);
vz = vz - mean(vz);

ax_old = 0*x;
ay_old = 0*y;
az_old = 0*z;
Ek = zeros(1,Nt);
Ep = zeros(1,Nt);

%% Verlet cell parameters
raw_cell_width = 2*G*D; % make cell size at least twice size of largest interaction distance (largest particle is D*G)
num_cells = (Lx/raw_cell_width);
num_cells = round(num_cells);
cell_width = Lx/num_cells;
jj_max = ceil(Lx/cell_width); % the maximum number of cells
cell_num = ceil(x/cell_width);

% create array (cell_list) where each particle is binned into a cell
for jj = 1:jj_max
    cell_list{jj} = find(cell_num == jj);
end
%% Setup Plotting
if plotit
    figure(1), clf;
    h=zeros(1,2*N);
    for np = 1:N
        h(np) = rectangle('Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
        h(np+N)=rectangle('Position',[Lx Ly Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
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
            for np = 1:N
                set(h(np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)]);
            end
            Np=N;
            ii=find(y<Dn/2);
            for nn=1:length(ii)
                np=ii(nn);
                set(h(nn+Np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np)+Ly Dn(np) Dn(np)]);
            end
            Np=Np+length(ii);
            %Top wall
            ii=find(y>Ly-Dn/2);
            for nn=1:length(ii)
                np=ii(nn);
                set(h(nn+Np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np)-Ly Dn(np) Dn(np)]);
            end
            Np=Np+length(ii);
            figure(1), ylim([0, Ly]);
            title(num2str(Ly));

        end
        figure(2), semilogy(nt, Ek(nt-1), 'ro');
        hold on, semilogy(nt, Ep(nt-1), 'bs');
        hold on, plot(nt,P,'kx')

%         figure(3), plot(nt,Ly,'ro',nt,Ly_min,'bx',nt,Ly_max,'kp'), hold on
        drawnow;

    elseif (plot_KE && mod(nt,Nplotskip) == 0)
        figure(3), plot3(x,y,z,'k.')
        axis equal
        drawnow
        % pause
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% First step in Verlet integration %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % move the particles according to last-loop's velocity and acceleration
    x  =  x+vx*dt+ax_old.*dt.^2/2;
    y  =  y+vy*dt+ay_old.*dt.^2/2;
    z  =  z+vz*dt+az_old.*dt.^2/2;

    % wrap around periodic boundaries
    y = mod(y,Ly); % !! this is not present in the 2d code 
    z = mod(z,Lz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Re-assign particles to cells %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if flag2 == true || mod(nt, Ncellupdate) == 0
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
    Fz = zeros(1,N);
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
                    dz = z(mm)-z(nn);
                    dz = dz - round(dz/Lz)*Lz;
                
                    dx = x(mm)-x(nn);
                    dx = dx - round(dx/Lx)*Lx;
                    dnm = dx.^2+dy.^2+dz.^2;
                    if(dnm < Dnm^2)
                        dnm = sqrt(dnm);

                        F = -K*(Dnm/dnm-1);

                        % dissipation force = B * m_red * v_N, v_N is normal component of velocity diff
                        m_red = M*M/(M+M);
                        v_dot_r=(vx(nn)-vx(mm))*dx + (vy(nn)-vy(mm))*dy + (vz(nn)-vz(mm))*dz;
                        Fdiss = Bv * m_red * v_dot_r;

                        Fx(nn) = Fx(nn)+F.*dx-Fdiss.*dx/dnm;  % particle-particle Force Law
                        %                         Fx(mm) = Fx(mm)-F.*dx+Fdiss.*dx/dnm;
                        Fy(nn) = Fy(nn)+F.*dy-Fdiss.*dy/dnm;
                        %                         Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
                        Fz(nn) = Fz(nn)+F.*dz-Fdiss.*dz/dnm;
                        Zn(nn) = Zn(nn) + 1;
                        %                         Zn(mm) = Zn(mm) + 1;
                        Ep(nt) = Ep(nt) + 0.5*K*(Dnm-dnm)^2;
                    end
                end
            end
        end

    end

    Fx = Fx - B.*vx;
    Fy = Fy - B.*vy;
    Fz = Fz - B.*vz;

    LW_contacts = x<Dn/2;
    % Fx = Fx-K*(x-Dn/2).*(LW_contacts);  % Left wall
    % Fy = Fy-K*(y-D/2).*(y<D/2);  % Bottom wall

    RW_contacts = x>Lx-Dn/2;
    % Fx = Fx-K*(x-(Lx-Dn/2)).*(RW_contacts);  % Right wall
    % Fy = Fy-K*(y-(Ly-D/2)).*(y>Ly-D/2);  % Top wall
    
    y=mod(y,Ly); %periodic boundaries for top and bottom
    x=mod(x,Lx); %periodic boundaries for left and right
    z=mod(z,Lz); %periodic boundaries for left and rightj

    Ek(nt) = 1/2*M*sum((vx).^2+(vy).^2+(vz).^2);
    Ek(nt) = Ek(nt)/N;
    Ep(nt) = Ep(nt)/N;


    ax = Fx./M;
    ay = Fy./M-g;
    az = Fz./M;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Second step in Verlet integration %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vx = vx+(ax_old+ax).*dt/2;
    vy = vy+(ay_old+ay).*dt/2;
    vz = vz+(az_old+az).*dt/2;
    
%     no_cont_list = (Zn == 0 & ~LW_contacts & ~RW_contacts);
%     vx(no_cont_list) = 0;
%     vy(no_cont_list) = 0;
%     ax(no_cont_list) = 0;
%     ay(no_cont_list) = 0;
% 

    ax_old = ax;
    ay_old = ay;
    az_old = az;

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

% save the tile
mass = M;
diameter_average = mean(Dn);
save(filename_tile, 'x', 'y', 'z', 'Dn', 'Lx', 'Ly', 'Lz', 'K', 'P_target', 'P', 'N', 'W_factor', 'mass', "diameter_average");

N_repeated = x_mult;
x_repeated = [];
y_repeated = [];
z_repeated = [];
Dn_repeated = [];

for i = 0:N_repeated-1
    x_shifted = x + i * Lx;
    y_shifted = y + i * Ly;
    z_shifted = z + i * Lz;
    x_repeated = [x_repeated, x_shifted];
    y_repeated = [y_repeated, y];
    z_repeated = [z_repeated, z];
    Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
end

N = N * N_repeated;
x = x_repeated;
y = y_repeated;
z = z_repeated;
Dn = Dn_repeated; % Set Dn to the repeated diameters

% Time for y
N_repeated = y_mult;
x_repeated = [];
y_repeated = [];
z_repeated = [];
Dn_repeated = [];

for i = 0:N_repeated-1
    y_shifted = y + i * Ly;
    x_repeated = [x_repeated, x];
    y_repeated = [y_repeated, y_shifted];
    z_repeated = [z_repeated, z];
    Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
end

N = N * N_repeated;
x = x_repeated;
y = y_repeated;
z = z_repeated;
Dn = Dn_repeated; % Set Dn to the repeated diameters
W_factor = W_factor * N_repeated;

% Time for z
N_repeated = z_mult;
x_repeated = [];
y_repeated = [];
z_repeated = [];
Dn_repeated = [];

for i = 0:N_repeated-1
    z_shifted = z + i * Lz;
    x_repeated = [x_repeated, x];
    y_repeated = [y_repeated, y];
    z_repeated = [z_repeated, z_shifted];
    Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
end

N = N * N_repeated;
x = x_repeated;
y = y_repeated;
z = z_repeated;
Dn = Dn_repeated; % Set Dn to the repeated diameters
% W_factor = W_factor * N_repeated; % !! Commented out because keeping the "width" a square

% figure;
% hold on;
% axis equal;
% for np = 1:N
%     rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], 'Curvature', [1, 1], 'EdgeColor', 'b');
% end
% axis([0, N_repeated * Lx, 0, Ly]);
% hold off;
% pause
% figure;

% hold on;
% axis equal;
% axis([0, x_mult * Lx, 0, Ly*y_mult]);

% % Loop through each particle and plot its rectangle
% % for np = 1:N
% %     % Compute the position for the rectangle using the center coordinates (x, y)
% %     % and the diameter Dn (width and height).
% %     rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], ...
% %         'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5); % Optional LineWidth for better visibility
% % end
% % Update the figure display
% drawnow;
% hold off;

% disp(['number of excess contacts = ' num2str(sum(Zn)/2 + sum(LW_contacts) + sum(RW_contacts) - 2*N)])

Lx = Lx * x_mult; 
Ly = Ly * y_mult;
Lz = Lz * z_mult;
W_factor = W_factor * y_mult; % could be either z or y, but keeping them the same for now
mass = M;
diameter_average = mean(Dn);

if calc_eig == true
    positions = [x',y'];
    radii = Dn./2;
    [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
    Hessian = hess2d(positions, radii, K, Ly, Lx);
    [eigen_vectors, eigen_values ] =  eig(Hessian);
    save(filename, 'x', 'y', 'Dn', 'Lx', 'Ly','Lz', 'K', 'P_target', 'P', 'N', 'eigen_vectors', 'eigen_values');
else
    save(filename, 'x', 'y', 'z', 'Dn', 'Lx', 'Ly', 'Lz', 'K', 'P_target', 'P', 'N', 'W_factor', 'mass', "diameter_average");
end

