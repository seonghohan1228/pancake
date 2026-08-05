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

%% Data
load_exp = [2.001499678
3.001444845
4.000304322
4.99854693
6.002785508
7.000460596
7.99801231
9.001584669
10.0049843
];
e_ratio_exp = [0.448205568
0.535694636
0.595467244
0.639491864
0.676586931
0.706123401
0.732510274
0.752597513
0.768275316
];

load = linspace (2, 10, 9);
e_ratio = [0.29, 0.42, 0.51, 0.57, 0.61, 0.65, 0.68, 0.70, 0.72]; % [um]

% Plot the data
figure;

% y-axis unit in kPa
hold on;
box on;
scatter(load_exp, e_ratio_exp, 'Marker', '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Experimental$^{\mathrm{(1)}}$');
scatter(load, e_ratio, 'Marker', 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'Numerical');
xlabel('$W\ \mathrm{(kN)}$');
ylabel('$\varepsilon$');
legend;
xlim([1.5 10.5]);
ylim([0.15 0.8]);
hold off;
