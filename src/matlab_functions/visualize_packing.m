load("in/2d_poly_20by20/2D_N12000_P0.001_Width20_Seed1.mat")

%  Below as been adddded ot it's own function called visualizeTile.m
% x_mult = 3;
% y_mult = 3;
% N_original = N;
% Dn_original = Dn;
% x_original = x;
% y_original = y;

% N_repeated = x_mult;
% x_repeated = [];
% y_repeated = [];
% Dn_repeated = [];

% for i = 0:N_repeated-1
%     x_shifted = x + i * Lx;
%     y_shifted = y + i * Ly;
%     x_repeated = [x_repeated, x_shifted];
%     y_repeated = [y_repeated, y];
%     Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
% end

% N = N * N_repeated;
% x = x_repeated;
% y = y_repeated;
% Dn = Dn_repeated; % Set Dn to the repeated diameters

% % Time for y
% N_repeated = y_mult;
% x_repeated = [];
% y_repeated = [];
% Dn_repeated = [];

% for i = 0:N_repeated-1
%     y_shifted = y + i * Ly;
%     x_repeated = [x_repeated, x];
%     y_repeated = [y_repeated, y_shifted];
%     Dn_repeated = [Dn_repeated, Dn]; % Append Dn for each repetition
% end

% N = N * N_repeated;
% x = x_repeated;
% y = y_repeated;
% Dn = Dn_repeated; % Set Dn to the repeated diameters
% % figure;
% % hold on;
% % axis equal;
% % for np = 1:N
% %     rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], 'Curvature', [1, 1], 'EdgeColor', 'b');
% % end
% % axis([0, N_repeated * Lx, 0, Ly]);
% % hold off;
% % pause
% figure;
% hold on;
% axis equal;
% axis([Lx, 2* Lx, Ly, Ly*2]);

% % Loop through each particle and plot its rectangle
% for np = 1:N
%     % Compute the position for the rectangle using the center coordinates (x, y)
%     % and the diameter Dn (width and height).
%     rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], ...
%         'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5); % Optional LineWidth for better visibility
    

%     % Wrap the "cell" index around the total number of particles
%     label = mod(np-1, N_original) + 1;  % This will cycle labels from 1 to N
%     text(x(np), y(np), num2str(label), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
% end
% % Update the figure display
% drawnow;
% hold off;

figure
hold on

% Loop through each particle and plot its rectangle
for np = 1:N_original
    rectangle('Position', [x_original(np) - Dn_original(np)/2, y_original(np) - Dn_original(np)/2, Dn_original(np), Dn_original(np)], ...
        'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5); 
end

box on
grid on
% xlim([1, 4])
% ylim([1,3])

% Get the start and end points for the arrows
x_start = x_original(3);
y_start = y_original(3);
% annotate only takes normalized values relative to figure limits. so need to convert data coordinates to normalized figure coordinates
ax = gca;
x_lim = ax.XLim; % X-axis limits [xmin, xmax]
y_lim = ax.YLim; 


% ====== X-direction arrow (Blue, pointing right) ======
x_end_x = x_start + Dn_original(3); % Arrow points in x-direction
y_end_x = y_start; 
x_start_x = (x_start - x_lim(1)) / (x_lim(2) - x_lim(1)); % Normalize the x-direction arrow
y_start_x = (y_start - y_lim(1)) / (y_lim(2) - y_lim(1));
x_end_x = (x_end_x - x_lim(1)) / (x_lim(2) - x_lim(1));
y_end_x = (y_end_x - y_lim(1)) / (y_lim(2) - y_lim(1));

% Draw annotation arrow
annotation('textarrow', [x_start_x x_end_x], [y_start_x y_end_x], 'Color', 'blue', 'Interpreter', 'Latex', 'FontSize', 12, 'LineWidth', 3);
% text(3, 2.3, '$A_x$', 'Interpreter', 'Latex', 'FontSize', 14, 'Color', 'blue');
annotation('textbox', [x_end_x, y_end_x, 0.05, 0.05], 'String', '$A_x$', 'Interpreter', 'Latex', 'FontSize', 14, 'EdgeColor', 'none', 'Color', 'blue');
% ====== y-direction arrow (red, pointing up) ======
x_end_y = x_start ;
y_end_y = y_start+ Dn_original(3); % Arrow points in y-direction 
x_start_y = (x_start - x_lim(1)) / (x_lim(2) - x_lim(1)); % Normalize the x-direction arrow
y_start_y = (y_start - y_lim(1)) / (y_lim(2) - y_lim(1));
x_end_y = (x_end_y - x_lim(1)) / (x_lim(2) - x_lim(1));
y_end_y = (y_end_y - y_lim(1)) / (y_lim(2) - y_lim(1));

% Draw annotation arrow
annotation('textarrow', [x_start_y x_end_y], [y_start_y y_end_y], 'Color', 'red', 'Interpreter', 'Latex', 'FontSize', 12, 'LineWidth', 3);
annotation('textbox', [x_end_y, y_end_y, 0.05, 0.05], 'String', '$A_y$', 'Interpreter', 'Latex', 'FontSize', 14, 'EdgeColor', 'none', 'Color', 'red');


% -----------------Below moved to own function visualizeFullPacking.m--------------------
% %  Full Channel 
% load("in/2d_5wide_1000long/2D_N5000_P0.001_Width5_Seed1.mat")
% load("in/2d_poly_20by20/2D_N12000_P0.001_Width20_Seed1.mat")
% figure;
% hold on;
% axis equal;
% % axis([0, Lx/60, 0, Ly]);

% % Loop through each particle and plot its rectangle
% for np = 1:N
%     % Compute the position for the rectangle using the center coordinates (x, y)
%     % and the diameter Dn (width and height).
%     rectangle('Position', [x(np) - Dn(np)/2, y(np) - Dn(np)/2, Dn(np), Dn(np)], ...
%         'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5); % Optional LineWidth for better visibility
% end
% % Update the figure display
% drawnow;
% hold off;