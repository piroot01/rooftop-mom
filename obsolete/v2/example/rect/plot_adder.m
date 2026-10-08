% open the existing figure
figHandle = openfig('ka_10p0_rect_2to1_CMA_10log10_ABS_sorted.fig');

% get the current axes from that figure and 'hold' it
ax = gca; 
hold(ax, 'on');

% plot the additional traces
load benchmark_results_3.mat

plot(ax, 1:length(lambda), 10 * log10(abs(lambda)), '--x', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'FFT accelerated version 15x30');

load benchmark_results_22x44.mat;

plot(ax, 1:length(lambda), 10 * log10(abs(lambda)), '--d', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'FFT accelerated version 22x44');

load benchmark_results_29x58.mat;

plot(ax, 1:length(lambda), 10 * log10(abs(lambda)), '--d', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'FFT accelerated version 29x58');

% update the legend to include the new entry
legend(ax, 'show');

% save the figure
savefig(figHandle, 'ka_10p0_rect_2to1_CHA_10log10_ABS_sorted_updated.fig');