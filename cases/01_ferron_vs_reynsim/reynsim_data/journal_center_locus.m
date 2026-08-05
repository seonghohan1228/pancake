%% Pre-processing
% Read the CSV file
data = readtable('journal_center_locus.csv', 'VariableNamingRule', 'preserve');

% Extract the columns from the table
t = data{:, "Time"};
x = data{:, "avg(x)"}*1e6;
y = data{:, "avg(y)"}*1e6;

% Create a folder to store the frames (if it doesn't exist)
outputFolder = 'frames';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
% Set the pause interval between frames (in seconds)
frameInterval = 0.1;

% Determine limits for the axes for consistency
xLimits = [-152*1.1, 152*1.1];
yLimits = [-152*1.1, 152*1.1];

% Define the radius
radius = 152;

% Create a vector of angles from 0 to 2*pi
theta = linspace(0, 2*pi, 100);

% Calculate the x and y coordinates of the circle
x_circle = radius * cos(theta);
y_circle = radius * sin(theta);

%% Plot
% Loop through each time step to plot the evolving locus
nFrames = length(t);
figure;

for k = 1:nFrames
    % Plot the data up to the current frame (line only)
    plot(x(1:k), y(1:k), 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', 'r');
    hold on;

    plot(x_circle, y_circle, 'k-', 'LineWidth', 1);
    % Add a filled circular marker at the current position only
    plot(x(k), y(k), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    
    % Labeling and formatting
    xlabel('$x\ \mathrm{(\mu m)}$', 'Interpreter', 'latex');
    ylabel('$y\ \mathrm{(\mu m)}$', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex', 'FontSize', 14);
    xlim(xLimits);
    ylim(yLimits);
    xticks([-150, -100, -50, 0, 50, 100, 150]);
    yticks([-150, -100, -50, 0, 50, 100, 150]);
    box on;
    axis square;
    
    hold off;
    drawnow;
    
    % Save the current figure as a PNG image at 300 DPI
    print(gcf, fullfile(outputFolder, sprintf('frame_%04d.png', k)), '-dpng', '-r300');

    % Pause for a consistent time interval between frames
    pause(frameInterval);
end
