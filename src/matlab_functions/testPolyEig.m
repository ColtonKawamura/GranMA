% Define the system parameters
m = 1;  % Mass matrix (m = 1 for simplicity)
k = 10;  % Spring constant (stiffness matrix)
beta = 0.05;  % Damping 
N = 10; % number of particles

% Define the matrices K, D, and M
K = 2*k * eye(N);  % Stiffness matrix
% make the off diagnal terms -k
for i = 1:N-1
    K(i, i+1) = -k;
    K(i+1, i) = -k;
end
% K(1, end) = -k;
% K(end, 1) = -k

M = m * eye(N);  % Mass matrix

D = 2*beta * eye(N); 
for i = 1:N-1
    D(i, i+1) = -beta;
    D(i+1, i) = -beta;
end

% System matrix (K + i*omega*beta*D - omega^2*M) 
% system_matrix = K + 1i * beta * D - M;  %

% Solve for the eigenvalues using the polyeig function
% eigenvalues = polyeig(system_matrix)

K
D
M

% % Try this way
[eigenvectors, eigenvalues] = polyeig(K, 1i*D, -M)

% for ii = 1:length(eigenvalues)
%     figure;
%     plot(real(eigenvectors(:, ii)), '-o');
%     title(['Eigenmode ' num2str(ii)]);
%     xlabel('particle position');
%     ylabel('Mode displacement');
%     grid on;
%     % pause; 
% end