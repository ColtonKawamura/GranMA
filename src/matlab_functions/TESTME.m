function [H, Zn] = TESTME(positions, radii, k, L_y, Lx)
    % Computes the Hessian matrix for a 2D granular packing with Hooke's force law.
    %
    % positions: Nx2 matrix, where each row is [x, y] for a particle
    % radii: N-vector, radius of each particle
    % k: Hooke's spring constant
    % L_y: Length of the box in the y-direction (for periodic boundary conditions)
    %
    % Returns:
    % H: 2N x 2N Hessian matrix (second derivatives of the potential)

    % Number of particles
    N = size(positions, 1);
    
    % Initialize the Hessian matrix
    H = zeros(2 * N, 2 * N);
    Zn = zeros(N,1);
    % ID the wall particles
    left_wall_list = (positions(:,1)<radii');
    right_wall_list = (positions(:,1)>Lx-radii');
    Zn(left_wall_list|right_wall_list) = 1;

    %% Clean the rattlers
    % Count the number of contatns to clean the rattlers 
    for i = 1:N
        for j = i+1+N
            % Compute the distance vector and magnitude
            dx = positions(i, 1) - positions(j, 1);
            dy = positions(i, 2) - positions(j, 2);
            
            % Apply periodic boundary conditions in the y-direction
            dy = dy - round(dy / L_y) * L_y;
            
            r = sqrt(dx^2 + dy^2);  % Euclidean distance
            
            % Compute the overlap distance (if any)
            overlap = radii(i) + radii(j) - r;
            
            if overlap > 0  % Particles are in contact
                Zn(i) = Zn(i)+1;
                Zn(j) = Zn(j)+1; 
            end
        end
    end

    % Remove particles with 2 or fewer contacts
    to_keep = Zn > 2; 

    positions = positions(to_keep, :);
    radii = radii(to_keep);
    Zn = Zn(to_keep);
    N = size(positions, 1);  % Update the number of particles     

    % do it one more time to get rid of "hidden" rattlers
    for i = 1:N
        for j = i+1+N
            % Compute the distance vector and magnitude
            dx = positions(i, 1) - positions(j, 1);
            dy = positions(i, 2) - positions(j, 2);
            
            % Apply periodic boundary conditions in the y-direction
            dy = dy - round(dy / L_y) * L_y;
            
            r = sqrt(dx^2 + dy^2);  % Euclidean distance
            
            % Compute the overlap distance (if any)
            overlap = radii(i) + radii(j) - r;
            
            if overlap > 0  % Particles are in contact
                Zn(i) = Zn(i)+1;
                Zn(j) = Zn(j)+1; 
            end
        end
    end
    % Remove particles with 2 or fewer contacts
    to_keep = Zn > 2; 

    positions = positions(to_keep, :);
    radii = radii(to_keep);
    Zn = Zn(to_keep);
    N = size(positions, 1);  % Update the number of particles     

    % Loop over all pairs of particles
    for i = 1:N
        for j = i+1:N
            % Compute the distance vector and magnitude
            dx = positions(i, 1) - positions(j, 1);
            dy = positions(i, 2) - positions(j, 2);
            
            % Apply periodic boundary conditions in the y-direction
            dy = dy - round(dy / L_y) * L_y;
            
            r = sqrt(dx^2 + dy^2);  % Euclidean distance
            
            % Compute the overlap distance (if any)
            overlap = radii(i) + radii(j) - r;
            
            if overlap > 0  % Particles are in contact
                Zn(i) = Zn(i)+1;
                Zn(j) = Zn(j)+1;

                % Compute elements of the Hessian for this pair
                Kxx = k * (1 - (dx^2 / r^2));
                Kyy = k * (1 - (dy^2 / r^2));
                Kxy = k * (-dx * dy / r^2);
                
                % Indices in the global Hessian matrix
                idx_i = 2 * i - 1;  % x-index for particle i
                idy_i = 2 * i;      % y-index for particle i
                idx_j = 2 * j - 1;  % x-index for particle j
                idy_j = 2 * j;      % y-index for particle j
                
                % Adding to zero, not a big deal because they are never double-counted.
                H(idx_i, idx_i) = H(idx_i, idx_i) + Kxx; % sub-matrix ii, position [1,1]
                H(idy_i, idy_i) = H(idy_i, idy_i) + Kyy; % sub-matrix ii, position [2,1]
                H(idx_i, idy_i) = H(idx_i, idy_i) + Kxy; % sub-matrix ii, position [1,2]
                H(idy_i, idx_i) = H(idy_i, idx_i) + Kxy; % sub-matrix ii, position [2,2]
                
                H(idx_j, idx_j) = H(idx_j, idx_j) + Kxx; % sub-matrix jj, position [1,1]
                H(idy_j, idy_j) = H(idy_j, idy_j) + Kyy; % sub-matrix jj, position [2,1]
                H(idx_j, idy_j) = H(idx_j, idy_j) + Kxy; % sub-matrix jj, position [1,1]
                H(idy_j, idx_j) = H(idy_j, idx_j) + Kxy; % sub-matrix jj, position [2,2]
                
                % Cross terms between particles i and j
                H(idx_i, idx_j) = H(idx_i, idx_j) - Kxx; % sub-matrix ij, position [1,1]
                H(idy_i, idy_j) = H(idy_i, idy_j) - Kyy; % sub-matrix ij, position [2,1]
                H(idx_i, idy_j) = H(idx_i, idy_j) - Kxy; % sub-matrix ij, position [1,2]
                H(idy_i, idx_j) = H(idy_i, idx_j) - Kxy; % sub-matrix ij, position [2,2]
                
                H(idx_j, idx_i) = H(idx_j, idx_i) - Kxx; % sub-matrix ji, position [1,1]
                H(idy_j, idy_i) = H(idy_j, idy_i) - Kyy; % sub-matrix ji, position [2,1]
                H(idx_j, idy_i) = H(idx_j, idy_i) - Kxy; % sub-matrix ji, position [1,2]
                H(idy_j, idx_i) = H(idy_j, idx_i) - Kxy; % sub-matrix ji, position [2,2]
            end
        end
    end
    for i = 1:N % go throuch each particle 
        if left_wall_list(i) || right_wall_list(i) % if the particle is either left or right wall
            idx_i = 2 * i - 1;  % x index for particle i
            idy_i = 2 * i;  % y index for particle i
            H(idx_i, idx_i) = H(idx_i, idx_i) + k;  % Add K to the x-coordinate diagonal Row ix_n = n and column ix_n = n correspond to the second derivative of the potential energy with respect to the x-coordinate of particle n
            H(idy_i, idy_i) = H(idy_i, idy_i) + k;  % same for y-coordinate diagonal
        end
    end
end

