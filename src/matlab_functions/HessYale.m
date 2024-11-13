function [Hessian, eigen_values, eigen_vectors] = HessYale(x, y, Dn, N, L, K, Lx)
% L = Ly (periodic in y)
% x = xy(1:N);
% y = xy(N+1:end);

% ID the wall particles
left_wall_list = (x<Dn/2);
right_wall_list = (x>Lx-Dn/2);
d2Ddr1dr2 = zeros(2*N, 2*N);

rad = Dn/2;

for n = 1:N-1 % start with the first particle (not sorted)
    ix_n = n; % index for the x-coord of particle n 
    iy_n = n + N; % index for y-coord
    for m = n+1:N % These are the neighbors of particle n
        ix_m = m; % index for partcile m's x-coord
        iy_m = m + N; % index y-coord
        Dnm = rad(n) + rad(m); % distance between n and m
        
        dx = x(m) - x(n);
        % dx = dx - round(dx / L) * L; % for periodic boundary in x
        if abs(dx) < Dnm
            dy = y(m) - y(n);
            dy = dy - round(dy / L) * L; % for periddib boundarty in y
            d = sqrt(dx^2 + dy^2);
            if d < Dnm
                dd = 1 - Dnm / d;
                Dnm_sq = Dnm^2;
                dx = dx / d; %
                dy = dy / d;
                dx_sq = dx^2;
                dy_sq = dy^2;
                % d2Ddx2 = (dx_sq + dd * dy_sq) / Dnm_sq / 2;
                % d2Ddy2 = (dy_sq + dd * dx_sq) / Dnm_sq / 2;
                % d2Ddxdy = dx * dy * (1 - dd) / Dnm_sq;

                d2Ddx2 = (dx_sq + dd * dy_sq) / 1/ 2;
                d2Ddy2 = (dy_sq + dd * dx_sq) / 1 / 2;
                d2Ddxdy = dx * dy * (1 - dd) / 1;

                d2Ddr1dr2(ix_n, ix_n) = d2Ddr1dr2(ix_n, ix_n) + d2Ddx2; % for n = 1: M(1, 1)
                d2Ddr1dr2(ix_n, iy_n) = d2Ddr1dr2(ix_n, iy_n) + d2Ddxdy; % M(1, 5001)
                d2Ddr1dr2(iy_n, iy_n) = d2Ddr1dr2(iy_n, iy_n) + d2Ddy2; %M(5001, 5001)

                d2Ddr1dr2(ix_m, ix_m) = d2Ddr1dr2(ix_m, ix_m) + d2Ddx2; 
                d2Ddr1dr2(ix_m, iy_m) = d2Ddr1dr2(ix_m, iy_m) + d2Ddxdy;
                d2Ddr1dr2(iy_m, iy_m) = d2Ddr1dr2(iy_m, iy_m) + d2Ddy2;

                d2Ddr1dr2(ix_n, ix_m) = d2Ddr1dr2(ix_n, ix_m) - 2 * d2Ddx2;
                d2Ddr1dr2(ix_n, iy_m) = d2Ddr1dr2(ix_n, iy_m) - d2Ddxdy;
                d2Ddr1dr2(ix_m, iy_n) = d2Ddr1dr2(ix_m, iy_n) - d2Ddxdy;
                d2Ddr1dr2(iy_n, iy_m) = d2Ddr1dr2(iy_n, iy_m) - 2 * d2Ddy2;

            end
        end
    end                    
end
% K = .01;
%% Adds K to the diagnals of the wall particles for non-perioldic sims
for i = 1:N % go throuch each particle 
    if left_wall_list(i) || right_wall_list(i) % if the particle is either left or right wall
        ix_i = i;  % x index for particle i
        iy_i = i + N;  % y index for particle i
        d2Ddr1dr2(ix_i, ix_i) = d2Ddr1dr2(ix_i, ix_i) + K/2;  % Add K to the x-coordinate diagonal Row ix_n = n and column ix_n = n correspond to the second derivative of the potential energy with respect to the x-coordinate of particle n
        d2Ddr1dr2(iy_i, iy_i) = d2Ddr1dr2(iy_i, iy_i) + K/2;  % same for y-coordinate diagonal
    end
end


Hessian = d2Ddr1dr2 + d2Ddr1dr2';

[eigen_vectors, eigen_values] = eig(Hessian); % [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.