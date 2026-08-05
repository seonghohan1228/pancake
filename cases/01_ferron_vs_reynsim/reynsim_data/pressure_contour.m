clc; clear; close all;

%% Apply figure styles
set(groot, ...
    'defaultAxesFontSize', 26, ...
    'defaultAxesFontName', 'Helvetica', ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesBox', 'on', ...
    'defaultLineMarkerFaceColor', 'w', ...
    'defaultLineLineWidth', 2, ...
    'defaultLineMarkerSize', 20, ... 
    'defaultFigureUnits', 'pixels', ...
    'defaultFigurePosition', [100 200 1000 700]);

%% Pressure contour
pressure_data = readtable('data/case_11_pressure.csv', 'VariableNamingRule', 'preserve');
p = pressure_data{:, "p"}/1e6;
z = pressure_data{:, "z"};
theta = pressure_data{:, "theta"}*180/pi;

% Define a grid for theta and z. Adjust nPts as needed.
nPtsTheta = 200;
nPtsZ = 50;
thetaGrid = linspace(min(theta), max(theta), nPtsTheta);
zGrid = linspace(min(z), max(z), nPtsZ);
[thetaMesh, zMesh] = meshgrid(thetaGrid, zGrid);

% Interpolate p onto the (theta, z) grid
pGrid = griddata(theta, z, p, thetaMesh, zMesh);

% Plot a filled contour
figure;
contourf(thetaGrid, zGrid, pGrid, 200, 'LineStyle', 'none');
colormap turbo;
colorbar;

xlabel('$\theta\ (^{\circ})$');
ylabel('$z\ \mathrm{(m)}$');

xlim([0 360]);
ylim([0 0.08]);
xticks([0, 90, 180, 270, 360]);
yticks([0, 0.02, 0.04, 0.06, 0.08]);

%% Pressure centerline
pressure_data = readtable('data/case_00_pressure_centerline.csv', 'VariableNamingRule', 'preserve');
p = pressure_data{:, "p"}/1e6;
theta = pressure_data{:, "theta"}*180/pi;

% Ferron 1983
pressure_exp_data = readtable('data/Ferron_pressure_centerline.csv', 'VariableNamingRule', 'preserve');
theta_exp = pressure_exp_data{:, "Var1"};
p_exp = pressure_exp_data{:, "Var2"};
% Shift to match angle start position to positive x-axis
theta_exp = mod((theta_exp + 90), 360);

% Plot a filled contour
figure; hold on;
plot(theta, p);
scatter(theta_exp, p_exp);
legend;

xlabel('$\theta\ (^{\circ})$');
ylabel('$p\ \mathrm{(MPa)}$');

xlim([0 360]);
xticks([0, 90, 180, 270, 360]);