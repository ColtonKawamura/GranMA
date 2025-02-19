function sim3d(K, M, Bv, w_D, N, P, W, seed, in_path, out_path)
    % sim3d(100, 1, 1, 10, 1060, 0.1, 5, 1, 'in/3D/3d_junkyard', 'out/simulation_3d/initial_test/')
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

% % % Manual variables for troubleshooting
% K = 100;
% M = 1;
% Bv = 1;
% w_D = 1.28;
% Nt =1000;
% N = 5000;
% P=0.1;
% W = 5;
% seed = 5;


plotdebug = 1;
% close all
K_in = K;
% Create the packing name with the exact number format for P
packing_name = string(sprintf("3D_N%d_P%s_Width%d_Seed%d", N, num2str(P), W, seed));

% Create the filename
% filename_output = sprintf('%s_K%d_Bv%d_wD%d_M%d.mat', packing_name, K, Bv, w_D, M);
filename_output = string(sprintf("%s_K%d_Bv%d_wD%.2f_M%d.mat", packing_name, K, Bv, w_D, M));

% if exist(char("out/simulation_2d/K100_everything_smaller_dt/" + filename_output), 'file')
%     return
% end
filename = in_path + packing_name + ".mat";  % Concatenate path and filename

load(filename);

input_pressure = P_target;

if exist(filename_output)
    fprintf('*** Alert: Output for this already exists. ***\n')
    return
end

K = K_in;
% N = length(Dn);
flag = true;
Lx0 = Lx;
D = min(Dn);
Lz = W*D;

B=0;

    
dt = pi*sqrt(M/K)*0.05; %  was pi*sqrt(M/K)*0.05
c_0 = min(Dn).*sqrt(K/M);
Nt = round(.9.*(Lx ./ c_0)./(dt))

A = P_target/100;

dt = pi*sqrt(M/K)*0.05;
ax_old = 0*x;
ay_old = 0*y;
az_old = 0*z;
vx = 0*x;
vy = 0*y;
vz = 0*z;

Ek = zeros(1,Nt);
Ep = zeros(1,Nt);
g = 0;

x_all = zeros(length(x),Nt);
y_all = x_all;
z_all = x_all;

%% initial positions

x0 = x;
y0 = y;
z0 = z;

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
            dz = z(mm)-z(nn);
            dz = dz - round(dz/Lz)*Lz;

            dx = x(mm)-x(nn);
            dnm = dx.^2+dy.^2+dz.^2;
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
% P = 0;
for nt = 1:Nt


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% First step in Verlet integration %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_all(:,nt) = x;
    y_all(:,nt) = y;
    z_all(:,nt) = z;

    x  =  x+vx*dt+ax_old.*dt.^2/2;
    y  =  y+vy*dt+ay_old.*dt.^2/2;
    z  =  z+vz*dt+az_old.*dt.^2/2;

    x(left_wall_list) = x0(left_wall_list)+A*sin(w_D*dt*nt);
    y(left_wall_list) = y0(left_wall_list);
    z(left_wall_list) = z0(left_wall_list);
    x(right_wall_list) = x0(right_wall_list);
    y(right_wall_list) = y0(right_wall_list);
    z(right_wall_list) = z0(right_wall_list);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Interaction detector and Force Law %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fx = zeros(1,N);
    Fy = zeros(1,N);
    Fz = zeros(1,N);
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
            dz = z(mm)-z(nn);
            dz = dz - round(dz/Lz)*Lz;

            dx = x(mm)-x(nn);
            dnm = dx.^2+dy.^2+dz.^2;
            %                 if(dnm < Dnm^2)
            dnm = sqrt(dnm);

            F = -K*(Dnm/dnm-1);

            % dissipation force = B * m_red * v_N, v_N is normal component of velocity diff
            %                     m_red = M*M/(M+M);
            dvx = vx(nn)-vx(mm);
            dvy = vy(nn)-vy(mm);
            dvz = vz(nn)-vz(mm);

            Fx(nn) = Fx(nn)+F.*dx-Bv*dvx;  % particle-particle Force Law
            %                     Fx(mm) = Fx(mm)-F.*dx+Fdiss.*dx/dnm;
            Fy(nn) = Fy(nn)+F.*dy-Bv*dvy;
            %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
            Fz(nn) = Fz(nn)+F.*dz-Bv*dvz;
            %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
            Zn(nn) = Zn(nn) + 1;
            %                     Zn(mm) = Zn(mm) + 1;
            Ep(nt) = Ep(nt) + 0.5*K*(Dnm-dnm)^2;
            %                 end
            %             end
        end
    end
    %
    % Fx = Fx - B.*vx;
    % Fy = Fy - B.*vy;
    % Fz = Fz - B.*vz;

    Fx(left_wall_list) = 0;
    Fx(right_wall_list) = 0;


    %     Fx = Fx-K*(x-Dn/2).*(x<Dn/2);  % Left wall
    % Fy = Fy-K*(y-D/2).*(y<D/2);  % Bottom wall

    %     Fx = Fx-K*(x-(Lx-Dn/2)).*(x>Lx-Dn/2);  % Right wall
    % Fy = Fy-K*(y-(Ly-D/2)).*(y>Ly-D/2);  % Top wall

    %     y=mod(y,Ly); %periodic boundaries for top and bottom

    Ek(nt) = 1/2*M*sum((vx).^2+(vy).^2);
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

    ax_old = ax;
    ay_old = ay;
    az_old = az;
end

%%%%% Post Processing %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Damping Warning for Data Anaylsis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the damping ratio
zeta = Bv / (2 * sqrt(M * K)); % Damping ratio formula

%% Check for overdamping
if zeta > 1
    warning('The system is overdamped! Damping ratio (zeta) is %f, which is greater than 1.', zeta);
elseif zeta == 1
    fprintf('The system is critically damped.\n');
else
    fprintf('The system is underdamped.\n');
end


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
    process_gm_fft(driving_amplitude, time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, initial_distance_from_oscillation);

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

% process_gm_fft_freq_density(time_vector, index_particles, index_oscillating_wall, driving_amplitude, position_particles, initial_distance_from_oscillation)
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
    process_gm_fft(driving_amplitude, time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, initial_distance_from_oscillation);

attenuation_y = fitted_attenuation;
attenuation_fit_line_y = attenuation_fit_line;
wavenumber_y = wavenumber;
unwrapped_phase_vector_y = unwrapped_phase_vector;
wavespeed_y = driving_frequency*2*pi*sqrt(M/K)/(wavenumber*1);
initial_distance_from_oscillation_output_y_fft = initial_distance_from_oscillation_output;
amplitude_vector_y = amplitude_vector;
cleaned_particle_index_y = cleaned_particle_index;
% process_gm_fft_freq_density(time_vector, index_particles, index_oscillating_wall, driving_amplitude, position_particles, initial_distance_from_oscillation)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Z Direction Post Processing
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert simulation variables to meet function convention
time_vector = (1:Nt)*dt;
[~,index_particles] = sort(z0);
index_oscillating_wall = left_wall_list;
driving_frequency = w_D/6.2832;
driving_amplitude=A;
position_particles = z_all;
initial_distance_from_oscillation = x0;

% Perform fft fitting
[fitted_attenuation, wavenumber, attenuation_fit_line, initial_distance_from_oscillation_output, amplitude_vector, unwrapped_phase_vector, cleaned_particle_index] = ...
    process_gm_fft(driving_amplitude, time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, initial_distance_from_oscillation);

% Don't waste time if we didn't get enough data
if isempty(cleaned_particle_index)

    fprintf('Simulation P=%d, Omega=%d, Gamma=%d, Seed=%d did not detect attenuation\n', P, w_D, Bv, seed);

    return
end

attenuation_z = fitted_attenuation;
attenuation_fit_line_z = attenuation_fit_line;
wavenumber_z = wavenumber;
unwrapped_phase_vector_z = unwrapped_phase_vector;
wavespeed_z = driving_frequency*2*pi*sqrt(M/K)/(wavenumber*1);
initial_distance_from_oscillation_output_z_fft = initial_distance_from_oscillation_output;
amplitude_vector_z = amplitude_vector;
cleaned_particle_index_z = cleaned_particle_index;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let's make everything Dimensionless
diameter_average = 1; % This is only here until all the old packings are updated to have this as an output CAO 1JUN2024
mass_particle_average = 1;
attenuation_x_dimensionless = attenuation_x*diameter_average;
attenuation_y_dimensionless = attenuation_y*diameter_average;
attenuation_z_dimensionless = attenuation_z*diameter_average;
wavenumber_x_dimensionless = wavenumber_x*diameter_average;
wavenumber_y_dimensionless = wavenumber_y*diameter_average;
wavenumber_z_dimensionless = wavenumber_z*diameter_average;
driving_angular_frequency_dimensionless = w_D*sqrt(mass_particle_average/K);
gamma_dimensionless = Bv/sqrt(K*mass_particle_average);
pressure_dimensionless = P;
% Save the file
save_path = fullfile(out_path, filename_output);

save( save_path, 'gamma_dimensionless','time_vector', 'index_particles', 'attenuation_x_dimensionless', ...
    'attenuation_y_dimensionless', 'attenuation_z_dimensionless', 'wavenumber_x_dimensionless', 'wavenumber_y_dimensionless', 'wavenumber_z_dimensionless', 'wavespeed_x', ...
     'wavespeed_y', 'wavespeed_z', 'driving_angular_frequency_dimensionless', 'attenuation_fit_line_x', ...
     'initial_distance_from_oscillation_output_x_fft', 'initial_distance_from_oscillation_output_y_fft','initial_distance_from_oscillation_output_z_fft', ...
     'amplitude_vector_x', 'amplitude_vector_y', 'amplitude_vector_z', "pressure_dimensionless", "seed", 'index_particles', 'input_pressure', 'unwrapped_phase_vector_x', 'unwrapped_phase_vector_y', 'unwrapped_phase_vector_z' );
