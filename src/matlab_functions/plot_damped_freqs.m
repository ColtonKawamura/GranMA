% This script is prototyping used during the JAN2025 trip to NPS to plot the effects of changing damping and pressure on the eigenvvalues/modes.


figure
hold on

% pressures = ["0.1", "0.01", "0.001"];
pressures = ["0.001"];
damping_constants = [1, 0.01];
% titles = {'High Damping', 'Medium Damping', 'Low Damping'};

figure
hold on

for i = 1:length(damping_constants)
    damping_constant = damping_constants(i);
    [normalized_variable, marker_color] = normVarColor(damping_constants, damping_constant, true);
    for j = 1:length(pressures)
        pressure = pressures(j);
        
        % filename = sprintf("in/damped_eig_test/2D_N100_P%s_Width10_Seed1.mat", pressure);
        filename = sprintf("in/2d_eigen_mode_test/2D_N400_P%s_Width10_Seed1.mat", pressure);
        load(filename)
        positions = [x', y'];
        radii = Dn ./ 2;
        [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
        mass = 1;
        [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
        [eigen_vectors, eigen_values] = polyeig(matSpring, matDamp, matMass);
        
        % scatter(abs(imag(eigen_values)), -real(eigen_values), damping_constant*100, 'DisplayName', sprintf('P=%.3f, %s', pressure, titles{i})) % different sizes for different damping
        plot(abs(imag(eigen_values)), -real(eigen_values)./damping_constant, ".",'Color', marker_color, 'DisplayName', sprintf('P=%.3f, damping= %.4f', pressure, damping_constant./10)) % different colors for different damping

    end
end

xlabel('Frequency')
ylabel('Damping')
grid on
legend
hold off
set(gca, "xscale", "log")
set(gca, "yscale", "log")




% % ------ Previus code
% load("in/damped_eig_test/2D_N14_P0.1_Width3_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.1;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass)


% % High Damping
% load("in/damped_eig_test/2D_N100_P0.1_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.5;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass); 
% [eigen_vectors_high, eigen_values_high] = polyeig(matSpring, matDamp, matMass);

% load("in/damped_eig_test/2D_N100_P0.01_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.5;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
% [eigen_vectors_mid, eigen_values_mid] = polyeig(matSpring, matDamp, matMass);

% load("in/damped_eig_test/2D_N100_P0.001_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.5;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
% [eigen_vectors_low, eigen_values_low] = polyeig(matSpring, matDamp, matMass);

% figure
% scatter( abs(imag(eigen_values_high)), -real(eigen_values_high), 'r', "DisplayName", 'P=0.1')
% hold on
% scatter( abs(imag(eigen_values_mid)), -real(eigen_values_mid), 'g', "DisplayName", 'P=0.01')
% scatter( abs(imag(eigen_values_low)), -real(eigen_values_low), 'b', "DisplayName", 'P=0.001')
% xlabel('Frequency')
% ylabel('Damping')
% grid on
% title('High Damping')

% % Medium Damping
% load("in/damped_eig_test/2D_N100_P0.1_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.05;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass)
% [eigen_vectors_high, eigen_values_high] = polyeig(matSpring, matDamp, matMass)

% load("in/damped_eig_test/2D_N100_P0.01_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.05;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass)
% [eigen_vectors_mid, eigen_values_mid] = polyeig(matSpring, matDamp, matMass)

% load("in/damped_eig_test/2D_N100_P0.001_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.05;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
% [eigen_vectors_low, eigen_values_low] = polyeig(matSpring, matDamp, matMass);

% figure
% scatter( abs(imag(eigen_values_high)), -real(eigen_values_high), 'r', "DisplayName", 'P=0.1')
% hold on
% scatter( abs(imag(eigen_values_mid)), -real(eigen_values_mid), 'g', "DisplayName", 'P=0.01')
% scatter( abs(imag(eigen_values_low)), -real(eigen_values_low), 'b', "DisplayName", 'P=0.001')
% xlabel('Frequency')
% ylabel('Damping')
% grid on
% title('Medium Damping')

% % Low Damping
% load("in/damped_eig_test/2D_N100_P0.1_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.005;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
% [eigen_vectors_high, eigen_values_high] = polyeig(matSpring, matDamp, matMass);

% load("in/damped_eig_test/2D_N100_P0.01_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.005;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
% [eigen_vectors_mid, eigen_values_mid] = polyeig(matSpring, matDamp, matMass);

% load("in/damped_eig_test/2D_N100_P0.001_Width10_Seed1.mat")
% positions = [x',y'];
% radii = Dn./2;
% [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% mass = 1;
% damping_constant = 0.005;
% [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
% [eigen_vectors_low, eigen_values_low] = polyeig(matSpring, matDamp, matMass);

% figure
% scatter( abs(imag(eigen_values_high)), -real(eigen_values_high), 'r', "DisplayName", 'P=0.1')
% hold on
% scatter( abs(imag(eigen_values_mid)), -real(eigen_values_mid), 'g', "DisplayName", 'P=0.01')
% scatter( abs(imag(eigen_values_low)), -real(eigen_values_low), 'b', "DisplayName", 'P=0.001')
% xlabel('Frequency')
% ylabel('Damping')
% grid on
% title('Low Damping')


% % Automated 
% damping_constants = [0.2, 0.05, 0.005];
% titles = {'High Damping', 'Medium Damping', 'Low Damping'};

% figure
% hold on

% for i = 1:length(damping_constants)
%         damping_constant = damping_constants(i);
        
%         % High Damping
%         load("in/damped_eig_test/2D_N100_P0.1_Width10_Seed1.mat")
%         positions = [x', y'];
%         radii = Dn ./ 2;
%         [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
%         mass = 1;
%         [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
%         [eigen_vectors_high, eigen_values_high] = polyeig(matSpring, matDamp, matMass);

%         load("in/damped_eig_test/2D_N100_P0.01_Width10_Seed1.mat")
%         positions = [x', y'];
%         radii = Dn ./ 2;
%         [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
%         mass = 1;
%         [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
%         [eigen_vectors_mid, eigen_values_mid] = polyeig(matSpring, matDamp, matMass);

%         load("in/damped_eig_test/2D_N100_P0.001_Width10_Seed1.mat")
%         positions = [x', y'];
%         radii = Dn ./ 2;
%         [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
%         mass = 1;
%         [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
%         [eigen_vectors_low, eigen_values_low] = polyeig(matSpring, matDamp, matMass);

%         scatter(abs(imag(eigen_values_high)), -real(eigen_values_high), damping_constant*100, 'r', "DisplayName", ['P=0.1, ' titles{i}])
%         scatter(abs(imag(eigen_values_mid)), -real(eigen_values_mid), damping_constant*100, 'g', "DisplayName", ['P=0.01, ' titles{i}])
%         scatter(abs(imag(eigen_values_low)), -real(eigen_values_low), damping_constant*100, 'b', "DisplayName", ['P=0.001, ' titles{i}])
% end

% xlabel('Frequency')
% ylabel('Damping')
% grid on
% legend
% title('Damping vs Frequency for Different Damping Constants')
% hold off


