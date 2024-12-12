function [positions, radii] = cleanRats(positions, radii, K, Ly, Lx)

    % Number of particles
    N = size(positions, 1);
    
    % Initialize the Hessian matrix
    Zn = zeros(N,1);
    % ID the wall particles
    left_wall_list = (positions(:,1)<radii');
    right_wall_list = (positions(:,1)>Lx-radii');
    Zn(left_wall_list|right_wall_list) = 2;
    %% Clean the rattlers
    % Count the number of contatns to clean the rattlers 
    for i = 1:N
        for j = i+1:N
            % Compute the distance vector and magnitude
            dx = positions(i, 1) - positions(j, 1);
            dy = positions(i, 2) - positions(j, 2);
            
            % Apply periodic boundary conditions in the y-direction
            dy = dy - round(dy / Ly) * Ly;
            
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
    % Zn = Zn(to_keep);
    N = size(positions,1);  % Update the number of particles     

    %% Re-intialize
    Zn = zeros(N,1);
    Zn(left_wall_list|right_wall_list) = 2;

    % do it one more time to get rid of "hidden" rattlers
    for i = 1:N
        for j = i+1:N
            % Compute the distance vector and magnitude
            dx = positions(i, 1) - positions(j, 1);
            dy = positions(i, 2) - positions(j, 2);
            
            % Apply periodic boundary conditions in the y-direction
            dy = dy - round(dy / Ly) * Ly;
            
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
    N = size(positions,1);  % Update the number of particles     
    %% Re-intialize
    Zn = zeros(N,1);
end