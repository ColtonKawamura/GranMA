% testTileDims.m


% Packings

% 20x20 tile (length times width)
save_path = 'in/2d_tile_test/20by20/'
packBi2dRepXY(400,100,1,1.4, 1, .1, 20, 1, false, 50, 1, false, save_path)
packBi2dRepXY(400,100,1,1.4, 1, .01, 20, 1, false, 50, 1, false, save_path)
packBi2dRepXY(400,100,1,1.4, 1, .001, 20, 1, false, 50, 1, false, save_path)

% 40 x 20
save_path = 'in/2d_tile_test/40by20/'
packBi2dRepXY(800,100,1,1.4, 1, .1, 20, 1, false, 25, 1, false, save_path)
packBi2dRepXY(800,100,1,1.4, 1, .01, 20, 1, false, 25, 1, false, save_path)
packBi2dRepXY(800,100,1,1.4, 1, .001, 20, 1, false, 25, 1, false, save_path)

% 100 x 20
save_path = 'in/2d_tile_test/100by20/'
packBi2dReXY(2000,100,1,1.4, 1, .1, 20, 1, false, 10, 1, false, save_path)
packBi2dRepXY(2000,100,1,1.4, 1, .01, 20, 1, false, 10, 1, false, save_path)
packBi2dRepXY(2000,100,1,1.4, 1, .001, 20, 1, false, 10, 1, false, save_path)

save_path = 'in/2d_tiled_2000by40/'
seed = [1,2,3,4,5]
pressure = [.1, .05, .01, .005, .001, .0005, .0001]
packBi2dRepXY(2000,100,1,1.4, 1, pressure, 20, seed, false, 20, 2, false, save_path)
packBi2dRepXY(2000,100,1,1.4, 1, .01, 20, 1, false, 20, 2, false, save_path)
packBi2dRepXY(2000,100,1,1.4, 1, .001, 20, 1, false, 20, 2, false, save_path)


% simulations
% 20x 20
in_path = "in/2d_tile_test/20by20/";
out_path = "out/simulation_2d/tile_test/20by20/"
K = 100;
M = 1;
Bv_list = [.001, .01, .1, 1, 10, 100];
w_D_list = [.001, .01, .1, 1, 10, 100];
N = 20000;
W = 20;
P_list = [.1, .01, .001]
seed = 1;
% sim2d(K, M, Bv, w_D, N, P, W, seed, in_path, out_path)

% Loop through every combination of Bv and w_D
for i = 1:length(Bv_list)
    for j = 1:length(w_D_list)
        for k = 1:length(P_list)
            Bv = Bv_list(i);  % Get current value of Bv
            w_D = w_D_list(j);  % Get current value of w_D
            P = P_list(k);
            disp(['Running job with Bv = ', num2str(Bv), ' and w_D = ', num2str(w_D), ' P = ', num2str(P)]);
            
            sim2d(K, M, Bv, w_D, N, P, W, seed, in_path, out_path) 
        end
    end
end

% Debugging
N = 5000
Bv = .01
w_D = 10
K = 100
M = 1
P = .1
W = 5
seed = 1
