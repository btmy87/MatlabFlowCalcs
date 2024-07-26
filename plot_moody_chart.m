%% plot_moody_chart
% generate Moody chart for Darcy friction factor

% Generate grid
re = logspace(log10(500), 8, 500);
roughness = [0, 1e-6, 1e-5, 1e-4, 2e-4, 5e-4, ...
             0.001, 0.002, 0.005, 0.01, 0.015, 0.02, ...
             0.03, 0.04, 0.05];
[RE, ROUGH] = ndgrid(re, roughness);

% Calculate friction factors
F = friction_factor(RE, ROUGH, implicit=true);
F2 = friction_factor(RE, ROUGH, implicit=false);

% plot chart
figure(Name="MoodyChart", Units="inches", Position=[1,1,10, 7.5]);
tiledlayout(1, 1, Padding="compact");
ha = nexttile;hold on;
xlabel("Reynolds Number, Re=\rhoVd/\mu")
ylabel("Darcy Friction Factor");
title("Moody Diagram");

ha.LineStyleOrder = ["-", "--", "-."];
h = gobjects(1, length(roughness));
for i = length(roughness):-1:1
    h(i) = plot(re, F(:, i), ...
        DisplayName=sprintf("\\epsilon/D=%.2g", roughness(i)));
end

% uncomment to plot Swamee-Jain approximation partially transparent
% ha.ColorOrderIndex = 1;
% ha.LineStyleOrderIndex = 1;
% for i = length(roughness):-1:1
%     h2 = plot(re, F2(:, i), ...
%         DisplayName=sprintf("Approx: \\epsilon/D=%.2g", roughness(i)));
%     h2.Color = [h2.Color, 0.3];
% end

% approximation won't appear in legend, too much
legend(h, Location="southwest");

ha.XScale = "log";
ha.YScale = "log";
ha.XLim = re([1,end]);