%% Pre-processing
% Read the CSV file
data = readtable('pressure_distribution_4000rpm_6000N.csv', 'VariableNamingRule', 'preserve');
exp_data = readtable('Ferron_pressure_distribution_4000rpm_6000N.csv', 'VariableNamingRule', 'preserve');

% Extract the columns from the table
p_raw = data{:, "p"}/1e6;
theta_raw = data{:, "theta"}*180/pi;

p_exp = [0.0464785187609383
0.0878437027321786
0.225842275716395
0.517825653225198
0.816243061156515
1.12111941829587
1.35064635176521
1.72103419604938
1.83098252279446
1.54126177101131
0.941961656634088
0.247661164321948
0.0333944951080276
0.0256131890148039
0.0309458091171197];
theta_exp = [29.3170508451206
59.7304285781922
89.5058854018279
120.208121696704
142.267426745000
157.509350774268
172.652496737716
187.265968996047
202.080386340851
217.138110707519
232.411731431975
247.957664644646
270.781478032845
300.722993643218
331.390742338927];

% Shift theta by subtracting 90 degrees and wrap to [0, 360)
theta_shifted = mod(theta_exp + 90, 360);

% Sort the shifted angles and rearrange p accordingly
[theta_sorted, sortIdx] = sort(theta_shifted);
p_shifted = p_exp(sortIdx);

%% Plot
figure; hold on;
scatter(theta_sorted, p_shifted, 50, 'k^', 'DisplayName', 'Experiment$^{(1)}$');
plot(theta_raw, p_raw, 'b-','LineWidth', 2, 'DisplayName', 'Simulation');
hold off;

% Format
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
xticks([0, 90, 180, 270, 360]);
fontsize(14,"points");
ylabel('$p\ \mathrm{(MPa)}$', 'Interpreter', 'latex');
set(gca,"TickLabelInterpreter",'latex');

hl = legend('show');
set(hl, 'Interpreter','latex', 'FontSize', 20);
set(gca,"TickLabelInterpreter",'latex', 'FontSize', 20);

xlim([0 360]);
ylim([0 2]);

box on;


%% Pressure contour
pressure_data = readtable('pressure_data_4000rpm_6000N.csv', 'VariableNamingRule', 'preserve');
p = pressure_data{:, "p"}/1e6;
z = pressure_data{:, "z"};
theta = pressure_data{:, "theta"}*180/pi;

% Define a grid for theta and z. Adjust nPts as needed.
nPtsTheta = 283;
nPtsZ = 72;
thetaGrid = linspace(min(theta), max(theta), nPtsTheta);
zGrid = linspace(min(z), max(z), nPtsZ);
[thetaMesh, zMesh] = meshgrid(thetaGrid, zGrid);

% Interpolate p onto the (theta, z) grid
pGrid = griddata(theta, z, p, thetaMesh, zMesh);

% Define the positions along the colormap scale (from 0 to 1)
positions = [0, 0.5, 1];

% Define the corresponding RGB colors at each position:
%   At 0: [0.2314, 0.2980, 0.7529]
%   At 0.5: [0.8650, 0.8650, 0.8650]
%   At 1: [0.7059, 0.0157, 0.1490]
colors = [0.23137254902, 0.298039215686, 0.752941176471; 
          0.8650,         0.8650,         0.8650;
          0.705882352941, 0.0156862745098, 0.149019607843];

% Create a fine grid for interpolation (e.g., 256 colors)
nColors = 256;
newPositions = linspace(0, 1, nColors);

% Interpolate to generate the full colormap
customMap = interp1(positions, colors, newPositions);

% Optional: Define a NaN color (e.g., yellow [1 1 0])
nanColor = [1, 1, 0];

% Plot a filled contour
figure;
contourf(thetaMesh, zMesh, pGrid, 20, 'LineStyle', 'none');
colormap(customMap);
caxis([0 1.8]);
% Create the colorbar and set its title
cb = colorbar;
cb.Title.String = '$p\ \mathrm(kPa)$';

% Optionally, you can set additional properties such as the font size or interpreter:
cb.Title.FontSize = 12;
cb.Title.Interpreter = 'latex';
% Set the tick label interpreter to LaTeX
cb.TickLabelInterpreter = 'latex';
xlabel('$\theta\ (^{\circ})$', 'Interpreter', 'latex');
ylabel('$z\ \mathrm{(m)}$', 'Interpreter', 'latex');

xlim([0 360]);
ylim([0 0.08]);
xticks([0, 90, 180, 270, 360]);
yticks([0, 0.02, 0.04, 0.06, 0.08]);
fontsize(20,"points");
set(gca,"TickLabelInterpreter",'latex');