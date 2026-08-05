W_exp = [1973.871065
2973.904431
3989.421271
5009.676943
6000.911384
7008.550739
7982.026902
9023.5415
10008.16543
] / 1000;

W_sim = [2000
3000
4000
5000
6000
7000
8000
9000
10000
] / 1000;

e_exp = [0.493289032
0.595543222
0.664413147
0.713713817
0.754386962
0.78009293
0.808683151
0.825749527
0.848582213
];
e_sim = [0.421
0.548
0.624
0.676
0.713
0.742
0.765
0.784
0.799 % Fix
];

% Plot the data
figure;
hold on;
scatter(W_exp, e_exp, 'k', '^', "filled"); % Use '-o' for a line with markers
scatter(W_sim, e_sim, 'k'); % Use '-o' for a line with markers
hold off;

% Set axis labels
xlabel('Load [kN]', 'FontSize', 14);
ylabel('Eccentricity ratio', 'FontSize', 14);

legend(["Ferron et al., 1983", "Simulation"]);

% Set x-axis limits to 0 and 360 degrees
xlim([1 11]);
ylim([0.2 0.9]);

% Add grid for clarity
grid on;
%title('Min. film thickness vs. Rotational speed', 'FontSize', 16);
