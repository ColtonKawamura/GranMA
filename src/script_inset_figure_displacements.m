% Open the figure and extract axes
openfig("shear_displacement_40wide_highP.fig")
fig_1 = gcf;
axObjs = findall(fig_1, 'Type', 'axes'); % Get all axes

% Extract data from each subplot
ax_left = axObjs(2);  % Left subplot
ax_right = axObjs(1); % Right subplot

left_dataObjs = ax_left.Children;
right_dataObjs = ax_right.Children;

leftXLim  = xlim(ax_left);
leftYLim = ylim(ax_left);

% Extract scatter data
xData_left = left_dataObjs(1).XData;
yData_left = left_dataObjs(1).YData;

xData_right = right_dataObjs(1).XData;
yData_right = right_dataObjs(1).YData;

% clean the rattlers for the visuliation
yData_left(yData_left == 0) = NaN;
yData_right(yData_right == 0) = NaN;

close(fig_1)

% Create a new figure because matlab is weird
figure;
ax1 = subplot(1,2,1); % Left subplot
ax2 = subplot(1,2,2); % Right subplot

% Replot data
scatter(ax1, xData_left, yData_left)
scatter(ax2, xData_right, yData_right)
xlim(ax1, [0, 1000])
xlim(ax2, [0, 1000])

box(ax1, 'on')
box(ax2, 'on')

xlabel(ax1, '$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel(ax1, '$\Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
xlabel(ax2, '$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel(ax2, '$\Delta y$', 'Interpreter', 'latex', 'FontSize', 20)

xlim(ax1, leftXLim)
ylim(ax1, leftYLim)
xlim(ax2, leftXLim)
ylim(ax2, leftYLim)

grid(ax1, 'on')
grid(ax2, 'on')
%% **Inset for the Left Subplot**
ax_inset1 = axes('Position', [0.3, 0.68, 0.19, 0.28]); % Adjust position & size
scatter(ax_inset1, xData_left, yData_left)
xlim(ax_inset1, [1, 100])
ylim(ax_inset1, [min(yData_left), max(yData_left)]) % Adjust y-axis for zoom-in
grid(ax_inset1, 'on')
set(ax_inset1, 'Box', 'on', 'Layer', 'top', 'Color', 'w')

%% **Inset for the Right Subplot**
ax_inset2 = axes('Position', [0.74, 0.68, 0.19, 0.28]); % Adjust position & size
scatter(ax_inset2, xData_right, yData_right)
xlim(ax_inset2, [1, 100])
ylim(ax_inset2, [min(yData_right), max(yData_right)]) % Adjust y-axis for zoom-in
grid(ax_inset2, 'on')
set(ax_inset2, 'Box', 'on', 'Layer', 'top', 'Color', 'w')

drawnow % Force MATLAB to refresh the figure

saveas(gcf, 'figures/shear_displacement_40wide_highP.eps', 'epsc');  % Save as color EPS

%% Shear Low pressure

% Open the figure and extract axes
openfig("shear_displacement_40wide_lowP.fig")
fig_1 = gcf;
axObjs = findall(fig_1, 'Type', 'axes'); % Get all axes

% Extract data from each subplot
ax_left = axObjs(2);  % Left subplot
ax_right = axObjs(1); % Right subplot

left_dataObjs = ax_left.Children;
right_dataObjs = ax_right.Children;

leftXLim  = xlim(ax_left);
leftYLim = ylim(ax_left);

% Extract scatter data
xData_left = left_dataObjs(1).XData;
yData_left = left_dataObjs(1).YData;

xData_right = right_dataObjs(1).XData;
yData_right = right_dataObjs(1).YData;

% clean the rattlers for the visuliation
yData_left(yData_left == 0) = NaN;
yData_right(yData_right == 0) = NaN;

close(fig_1)

% Create a new figure because matlab is weird
figure;
ax1 = subplot(1,2,1); % Left subplot
ax2 = subplot(1,2,2); % Right subplot

% Replot data
scatter(ax1, xData_left, yData_left)
scatter(ax2, xData_right, yData_right)

box(ax1, 'on')
box(ax2, 'on')

xlim(ax1, [0, 1000])
xlim(ax2, [0, 1000])

xlabel(ax1, '$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel(ax1, '$\Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
xlabel(ax2, '$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel(ax2, '$\Delta y$', 'Interpreter', 'latex', 'FontSize', 20)

xlim(ax1, leftXLim)
ylim(ax1, leftYLim)
xlim(ax2, leftXLim)
ylim(ax2, leftYLim)

grid(ax1, 'on')
grid(ax2, 'on')
box on
%% **Inset for the Left Subplot**
ax_inset1 = axes('Position', [0.3, 0.68, 0.19, 0.28]); % Adjust position & size
scatter(ax_inset1, xData_left, yData_left)
xlim(ax_inset1, [1, 40])
ylim(ax_inset1, [min(yData_left), max(yData_left)]) % Adjust y-axis for zoom-in
grid(ax_inset1, 'on')
set(ax_inset1, 'Box', 'on', 'Layer', 'top', 'Color', 'w')

%% **Inset for the Right Subplot**
ax_inset2 = axes('Position', [0.74, 0.68, 0.19, 0.28]); % Adjust position & size
scatter(ax_inset2, xData_right, yData_right)
xlim(ax_inset2, [1, 40])
ylim(ax_inset2, [min(yData_right), max(yData_right)]) % Adjust y-axis for zoom-in
grid(ax_inset2, 'on')
set(ax_inset2, 'Box', 'on', 'Layer', 'top', 'Color', 'w')

drawnow % Force MATLAB to refresh the figure
saveas(gcf, 'figures/shear_displacement_40wide_lowP.eps', 'epsc');  % Save as color EPS


%% Compression High pressure

% Open the figure and extract axes
openfig("compression_displacement_40wide_highP.fig")
fig_1 = gcf;
axObjs = findall(fig_1, 'Type', 'axes'); % Get all axes

% Extract data from each subplot
ax_left = axObjs(2);  % Left subplot
ax_right = axObjs(1); % Right subplot

left_dataObjs = ax_left.Children;
right_dataObjs = ax_right.Children;

leftXLim  = xlim(ax_left);
leftYLim = ylim(ax_left);

% Extract scatter data
xData_left = left_dataObjs(1).XData;
yData_left = left_dataObjs(1).YData;

xData_right = right_dataObjs(1).XData;
yData_right = right_dataObjs(1).YData;

% clean the rattlers for the visuliation
yData_left(yData_left == 0) = NaN;
yData_right(yData_right == 0) = NaN;

close(fig_1)

% Create a new figure because matlab is weird
figure;
ax1 = subplot(1,2,1); % Left subplot
ax2 = subplot(1,2,2); % Right subplot

% Replot data
scatter(ax1, xData_left, yData_left)
scatter(ax2, xData_right, yData_right)

box(ax1, 'on')
box(ax2, 'on')

xlim(ax1, [0, 1000])
xlim(ax2, [0, 1000])

xlabel(ax1, '$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel(ax1, '$\Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
xlabel(ax2, '$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel(ax2, '$\Delta y$', 'Interpreter', 'latex', 'FontSize', 20)

xlim(ax1, leftXLim)
ylim(ax1, leftYLim)
xlim(ax2, leftXLim)
ylim(ax2, leftYLim)

grid(ax1, 'on')
grid(ax2, 'on')
box on
%% **Inset for the Left Subplot**
ax_inset1 = axes('Position', [0.3, 0.68, 0.19, 0.28]); % Adjust position & size
scatter(ax_inset1, xData_left, yData_left)
xlim(ax_inset1, [1, 100])
ylim(ax_inset1, [min(yData_left), max(yData_left)]) % Adjust y-axis for zoom-in
grid(ax_inset1, 'on')
set(ax_inset1, 'Box', 'on', 'Layer', 'top', 'Color', 'w')

%% **Inset for the Right Subplot**
ax_inset2 = axes('Position', [0.74, 0.68, 0.19, 0.28]); % Adjust position & size
scatter(ax_inset2, xData_right, yData_right)
xlim(ax_inset2, [1, 100])
ylim(ax_inset2, [min(yData_right), max(yData_right)]) % Adjust y-axis for zoom-in
grid(ax_inset2, 'on')
set(ax_inset2, 'Box', 'on', 'Layer', 'top', 'Color', 'w')

drawnow % Force MATLAB to refresh the figure

saveas(gcf, 'figures/compression_displacement_40wide_highP.eps', 'epsc');  % Save as color EPS
%% Compression Low pressure

% Open the figure and extract axes
openfig("compression_displacement_40wide_lowP.fig")
fig_1 = gcf;
axObjs = findall(fig_1, 'Type', 'axes'); % Get all axes

% Extract data from each subplot
ax_left = axObjs(2);  % Left subplot
ax_right = axObjs(1); % Right subplot

left_dataObjs = ax_left.Children;
right_dataObjs = ax_right.Children;

leftXLim  = xlim(ax_left);
leftYLim = ylim(ax_left);

% Extract scatter data
xData_left = left_dataObjs(1).XData;
yData_left = left_dataObjs(1).YData;

xData_right = right_dataObjs(1).XData;
yData_right = right_dataObjs(1).YData;

close(fig_1)
% Create a new figure because matlab is weird
figure;
ax1 = subplot(1,2,1); % Left subplot
ax2 = subplot(1,2,2); % Right subplot

% clean the rattlers for the visuliation
yData_left(yData_left == 0) = NaN;
yData_right(yData_right == 0) = NaN;

% Replot data
scatter(ax1, xData_left, yData_left)
scatter(ax2, xData_right, yData_right)

box(ax1, 'on')
box(ax2, 'on')

xlim(ax1, [0, 1000])
xlim(ax2, [0, 1000])

xlabel(ax1, '$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel(ax1, '$\Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
xlabel(ax2, '$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel(ax2, '$\Delta y$', 'Interpreter', 'latex', 'FontSize', 20)

xlim(ax1, leftXLim)
ylim(ax1, leftYLim)
xlim(ax2, leftXLim)
ylim(ax2, leftYLim)

grid(ax1, 'on')
grid(ax2, 'on')
box on
%% **Inset for the Left Subplot**
ax_inset1 = axes('Position', [0.3, 0.68, 0.19, 0.28]); % Adjust position & size
scatter(ax_inset1, xData_left, yData_left)
xlim(ax_inset1, [1, 100])
ylim(ax_inset1, [min(yData_left), max(yData_left)]) % Adjust y-axis for zoom-in
grid(ax_inset1, 'on')
set(ax_inset1, 'Box', 'on', 'Layer', 'top', 'Color', 'w')

%% **Inset for the Right Subplot**
ax_inset2 = axes('Position', [0.74, 0.68, 0.19, 0.28]); % Adjust position & size
scatter(ax_inset2, xData_right, yData_right)
xlim(ax_inset2, [1, 100])
ylim(ax_inset2, [min(yData_right), max(yData_right)]) % Adjust y-axis for zoom-in
grid(ax_inset2, 'on')
set(ax_inset2, 'Box', 'on', 'Layer', 'top', 'Color', 'w')

drawnow % Force MATLAB to refresh the figure
saveas(gcf, 'figures/compression_displacement_40wide_lowP.eps', 'epsc');  % Save as color EPS