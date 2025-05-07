function [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, k, Ly, Lx, damping_constant, mass)
    % Computes the Hessian matrix for a 2D granular packing with Hooke's force law.
    %
    % positions: Nx2 matrix, where each row is [x, y] for a particle
    % radii: N-vector, radius of each particle
    % k: Hooke's spring constant
    % Ly: Length of the box in the y-direction (for periodic boundary conditions)
    %
    % Returns:
    % Hessian: 2N x 2N Hessian matrix (second derivatives of the potential)
    % k =1;
    % damping_constant = 1;
    N = size(positions,1);
    Zn = zeros(N,1);
    left_wall_list = (positions(:,1)<radii);
    right_wall_list = (positions(:,1)>Lx-radii);
    Zn(left_wall_list|right_wall_list) = 2;
    matSpring = zeros(2 * N, 2 * N);
    matDamp = zeros(2 * N, 2 * N);

    % Loop over all pairs of particles
    for i = 1:N
        for j = i+1:N

            % distnaces 
            dx = positions(i, 1) - positions(j, 1);
            dy = positions(i, 2) - positions(j, 2);

            % Aapply periodic boundary conditions
            dy = dy - round(dy / Ly) * Ly;
            dx= dx - round(dx / Lx) * Lx;
            r = sqrt(dx^2 + dy^2); 
            
            overlap = radii(i) + radii(j) - r;
            
            if overlap > 0 

                % keep tally of contacts for each partcile
                Zn(i) = Zn(i)+1;
                Zn(j) = Zn(j)+1;

                Kxx = k * ((dx^2 / r^2));
                Kyy = k * ((dy^2 / r^2));
                Kxy = k * (-dx * dy / r^2); % same as Kyx
                
                % Damping Matrix
                Dxx = 1;
                Dyy = 1;
                Dxy = -1;

                % Indices in the global Hessian matrix
                idx_i = 2 * i - 1;  % x-index for particle i
                idy_i = 2 * i;      % y-index for particle i
                idx_j = 2 * j - 1;  % x-index for particle j
                idy_j = 2 * j;      % y-index for particle j
                
                % Adding to zero, not a big deal because they are never double-counted.
                matSpring(idx_i, idx_i) = matSpring(idx_i, idx_i) + Kxx; % sub-matrix ii, position [1,1]
                matSpring(idy_i, idy_i) = matSpring(idy_i, idy_i) + Kyy; % sub-matrix ii, position [2,1]
                matSpring(idx_i, idy_i) = matSpring(idx_i, idy_i) + Kxy; % sub-matrix ii, position [1,2]
                matSpring(idy_i, idx_i) = matSpring(idy_i, idx_i) + Kxy; % sub-matrix ii, position [2,2]
                
                matSpring(idx_j, idx_j) = matSpring(idx_j, idx_j) + Kxx; % sub-matrix jj, position [1,1]
                matSpring(idy_j, idy_j) = matSpring(idy_j, idy_j) + Kyy; % sub-matrix jj, position [2,1]
                matSpring(idx_j, idy_j) = matSpring(idx_j, idy_j) + Kxy; % sub-matrix jj, position [1,1]
                matSpring(idy_j, idx_j) = matSpring(idy_j, idx_j) + Kxy; % sub-matrix jj, position [2,2]
                
                % Cross terms between particles i and j
                matSpring(idx_i, idx_j) = matSpring(idx_i, idx_j) - Kxx; % sub-matrix ij, position [1,1]
                matSpring(idy_i, idy_j) = matSpring(idy_i, idy_j) - Kyy; % sub-matrix ij, position [2,1]
                matSpring(idx_i, idy_j) = matSpring(idx_i, idy_j) - Kxy; % sub-matrix ij, position [1,2]
                matSpring(idy_i, idx_j) = matSpring(idy_i, idx_j) - Kxy; % sub-matrix ij, position [2,2]
                
                matSpring(idx_j, idx_i) = matSpring(idx_j, idx_i) - Kxx; % sub-matrix ji, position [1,1]
                matSpring(idy_j, idy_i) = matSpring(idy_j, idy_i) - Kyy; % sub-matrix ji, position [2,1]
                matSpring(idx_j, idy_i) = matSpring(idx_j, idy_i) - Kxy; % sub-matrix ji, position [1,2]
                matSpring(idy_j, idx_i) = matSpring(idy_j, idx_i) - Kxy; % sub-matrix ji, position [2,2]

                % Damping Matrix
                matDamp(idx_i, idx_i) = matDamp(idx_i, idx_i) + Dxx; % sub-matrix ii, position [1,1]
                matDamp(idy_i, idy_i) = matDamp(idy_i, idy_i) + Dyy; % sub-matrix ii, position [2,1]
                %matDamp(idx_i, idy_i) = matDamp(idx_i, idy_i) + Dxy; % sub-matrix ii, position [1,2]
               % matDamp(idy_i, idx_i) = matDamp(idy_i, idx_i) + Dxy; % sub-matrix ii, position [2,2]

                matDamp(idx_j, idx_j) = matDamp(idx_j, idx_j) + Dxx; % sub-matrix jj, position [1,1]
                matDamp(idy_j, idy_j) = matDamp(idy_j, idy_j) + Dyy; % sub-matrix jj, position [2,1]
                %matDamp(idx_j, idy_j) = matDamp(idx_j, idy_j) + Dxy; % sub-matrix jj, position [1,1]
                %matDamp(idy_j, idx_j) = matDamp(idy_j, idx_j) + Dxy; % sub-matrix jj, position [2,2]

                % Cross terms between particles i and j
                matDamp(idx_i, idx_j) = matDamp(idx_i, idx_j) - Dxx; % sub-matrix ij, position [1,1]
                matDamp(idy_i, idy_j) = matDamp(idy_i, idy_j) - Dyy; % sub-matrix ij, position [2,1]
                %matDamp(idx_i, idy_j) = matDamp(idx_i, idy_j) - Dxy; % sub-matrix ij, position [1,2]
                %matDamp(idy_i, idx_j) = matDamp(idy_i, idx_j) - Dxy; % sub-matrix ij, position [2,2]

                matDamp(idx_j, idx_i) = matDamp(idx_j, idx_i) - Dxx; % sub-matrix ji, position [1,1]
                matDamp(idy_j, idy_i) = matDamp(idy_j, idy_i) - Dyy; % sub-matrix ji, position [2,1]
                %matDamp(idx_j, idy_i) = matDamp(idx_j, idy_i) - Dxy; % sub-matrix ji, position [1,2]
                %matDamp(idy_j, idx_i) = matDamp(idy_j, idx_i) - Dxy; % sub-matrix ji, position [2,2]

            end
        end
    end
    for i = 1:N % go throuch each particle 
        if left_wall_list(i) || right_wall_list(i) % if the particle is either left or right wall
            idx_i = 2 * i - 1;  % x index for particle i
            idy_i = 2 * i;  % y index for particle i
            matSpring(idx_i, idx_i) = matSpring(idx_i, idx_i) + k;  % Add K to the x-coordinate diagonal Row ix_n = n and column ix_n = n correspond to the second derivative of the potential energy with respect to the x-coordinate of particle n
            matSpring(idy_i, idy_i) = matSpring(idy_i, idy_i) + k;  % same for y-coordinate diagonal
        end
    end
    matMass = mass*eye(2*N);
    matDamp = damping_constant*matDamp;
end

