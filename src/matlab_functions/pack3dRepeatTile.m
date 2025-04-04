function pack3dRepeatTile(N, K, P_target, W_factor, seed, x_mult, y_mult, z_mult, calc_eig, in_path, save_path)
% pack3dRepeatTile(3780, 100, .1, 15, 1, 40, 1, 1, false, 'in/3D/tile_15by15by15/', 'in/3D/tiled_15by600/')

packing_name = string(sprintf("3D_tile_N%d_P%s_Width%d_Seed%d", N, num2str(P_target), W_factor, seed));

filename = in_path + packing_name + ".mat";  % Concatenate path and filename

load(filename);

N_new = N * x_mult * y_mult * z_mult;
W_new = W_factor * y_mult; % !! Currently, this only inclues square-width packings (y_mult = z_mult)


N_repeated = x_mult;
x_repeated = [];
y_repeated = [];
z_repeated = [];
Dn_repeated = [];
filename = [save_path, '3D_N' num2str(N_new) '_P' num2str(P_target) '_Width' num2str(W_new) '_Seed' num2str(seed) '.mat'];

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

Lx = Lx * x_mult; 
Ly = Ly * y_mult;
Lz = Lz * z_mult;
W_factor = W_factor * y_mult; % could be either z or y, but keeping them the same for now
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