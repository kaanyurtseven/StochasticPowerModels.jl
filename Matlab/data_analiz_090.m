clear all
close all
clc

%%

load GenerationData090

gen = table2array(GenerationData090);

gen(:,2) = gen(:,2) + abs(min(gen(:,2)));



gen_det = [0.988962 0.14739];

gen_exp = [mean(gen(:,1)) mean(gen(:,2))];


%%
confidenceLevel = 0.05;

% Sort the data in ascending order
sortedData = sort(gen(:,1));

% Calculate the VaR
VaR_index = ceil((1 - confidenceLevel) * length(sortedData));
VaR(:,1) = sortedData(VaR_index);

% Calculate the CVaR
CVaR(:,1) = mean(sortedData(VaR_index:end));


%%
sortedData = sort(gen(:,2));

% Calculate the VaR
VaR_index = ceil((1 - confidenceLevel) * length(sortedData));
VaR(:,2) = sortedData(VaR_index);

% Calculate the CVaR
CVaR(:,2) = mean(sortedData(VaR_index:end));

gen_cvar = CVaR;


%%

gen_re_det1 = gen(:,1) - gen_det(:,1);
gen_re_det1(gen_re_det1 < 0) = [];

gen_re_det2 = 2*(gen(:,2) - gen_det(:,2));
gen_re_det2(gen_re_det2 < 0) = [];

gen_re_exp1 = gen(:,1) - gen_exp(:,1);
gen_re_exp1(gen_re_exp1 < 0) = [];

gen_re_exp2 = 2*(gen(:,2) - gen_exp(:,2));
gen_re_exp2(gen_re_exp2 < 0) = [];

gen_re_cvar1 = gen(:,1) - gen_cvar(:,1);
gen_re_cvar1(gen_re_cvar1 < 0) = [];

gen_re_cvar2 = 2*(gen(:,2) - gen_cvar(:,2));
gen_re_cvar2(gen_re_cvar2 < 0) = [];

exp_cost_det = mean(gen_re_det1) + mean(gen_re_det2)
exp_cost_exp = mean(gen_re_exp1) + mean(gen_re_exp2)
exp_cost_cvar = mean(gen_re_cvar1) + mean(gen_re_cvar2)


%% Visualization

color1 = [244, 104, 96] ./ 255;
color2 = [188, 143, 143] ./ 255;
color3 = [135, 206, 235] ./ 255;
color4 = [171, 196, 142] ./ 255;

f = figure;
f.Position = [100 100 1600 600];

bin = 50;

font = 14;
linewidth = 1.5;


t = tiledlayout(1,2,'TileSpacing','loose');

%Tile 1
nexttile
histogram(gen(:,1), bin, 'Normalization','pdf', 'FaceColor', color4, 'EdgeColor', color4)
xlabel('\it Gen 1', 'FontSize', font)
xline(gen_det(1), '--', '\it Det', 'LineWidth', linewidth, 'FontSize', font, 'Color', color2)
xline(gen_exp(1), '--', '\it Exp', 'LineWidth', linewidth, 'FontSize', font, 'Color', color3)
xline(gen_cvar(1), '--', '$CVaR_{0.95}$', 'LineWidth', linewidth, 'FontSize', font+2, 'Color', color1, 'Interpreter', 'latex');
xlim([0 2.2])
ylim([0 4])

nexttile
histogram(gen(:,2), bin, 'Normalization','pdf', 'FaceColor', color4, 'EdgeColor', color4)
xlabel('\it Gen 2', 'FontSize', font)
xline(gen_det(2), '--', '\it Det', 'LineWidth', linewidth, 'FontSize', font, 'Color', color2)
xline(gen_exp(2), '--', '\it Exp', 'LineWidth', linewidth, 'FontSize', font, 'Color', color3)
xline(gen_cvar(2), '--', '$CVaR_{0.95}$', 'LineWidth', linewidth, 'FontSize', font+2, 'Color', color1, 'Interpreter', 'latex');
xlim([0 2.2])
ylim([0 4])

xlabel(t, [newline '\it Active Power [pu]'], 'FontSize', font)
ylabel(t, '\it Density', 'FontSize', font)


%% 
f = figure;
f.Position = [100 100 1600 600];
t = tiledlayout(2,3,'TileSpacing','loose');
x_max1 = 0.7;
x_max2 = 2.5;
y_max1 = 12;
y_max2 = 5;

%Tile 1
nexttile
histogram(gen_re_det1, bin, 'Normalization','pdf', 'FaceColor', color4, 'EdgeColor', color4)
xlim([0 x_max1])
ylim([0 y_max1])
xlabel('\it Gen 1 - Det', 'FontSize', font)
xline(mean(gen_re_det1), '--', ['\it', num2str(round(mean(gen_re_det1),3))], 'LineWidth', linewidth, 'FontSize', font, 'Color', color2)
% L = legend('Det', 'Exp');

%Tile 2
nexttile
histogram(gen_re_exp1, bin, 'Normalization','pdf', 'FaceColor', color4, 'EdgeColor', color4)
xlim([0 x_max1])
ylim([0 y_max1])
xlabel('\it Gen 1 - Exp', 'FontSize', font)
xline(mean(gen_re_exp1), '--', ['\it', num2str(round(mean(gen_re_exp1),3))], 'LineWidth', linewidth, 'FontSize', font, 'Color', color3)

%Tile 3
nexttile
histogram(gen_re_cvar1, bin, 'Normalization','pdf', 'FaceColor', color4, 'EdgeColor', color4)
xlim([0 x_max1])
ylim([0 y_max1])
xlabel('\it Gen 1 - CVaR', 'FontSize', font)
xline(mean(gen_re_cvar1), '--', ['\it', num2str(round(mean(gen_re_cvar1),3))], 'LineWidth', linewidth, 'FontSize', font, 'Color', color1)


%Tile 1
nexttile
histogram(gen_re_det2, bin, 'Normalization','pdf', 'FaceColor', color4, 'EdgeColor', color4)
xlim([0 x_max2])
ylim([0 y_max2])
xlabel('\it Gen 2 - Det', 'FontSize', font)
xline(mean(gen_re_det2), '--', ['\it', num2str(round(mean(gen_re_det2),3))], 'LineWidth', linewidth, 'FontSize', font, 'Color', color2)

%Tile 2
nexttile
histogram(gen_re_exp2, bin, 'Normalization','pdf', 'FaceColor', color4, 'EdgeColor', color4)
xlim([0 x_max2])
ylim([0 y_max2])
xlabel('\it Gen 2 - Exp', 'FontSize', font)
xline(mean(gen_re_exp2), '--', ['\it', num2str(round(mean(gen_re_exp2),3))], 'LineWidth', linewidth, 'FontSize', font, 'Color', color3)

%Tile 3
nexttile
histogram(gen_re_cvar2, bin, 'Normalization','pdf', 'FaceColor', color4, 'EdgeColor', color4)
xlim([0 x_max2])
ylim([0 y_max2])
xlabel('\it Gen 2 - CVaR', 'FontSize', font)
xline(mean(gen_re_cvar2), '--', ['\it', num2str(round(mean(gen_re_cvar2),3))], 'LineWidth', linewidth, 'FontSize', font, 'Color', color1)

xlabel(t, [newline '\it Redispatch Unit Cost'], 'FontSize', font)
ylabel(t, '\it Density', 'FontSize', font)










