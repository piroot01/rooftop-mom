% MATLAB Benchmark Script
%
% This script compares the performance of three different methods for
% computing a masked matrix-vector product y = A*v, where parts of the
% matrix and vector are effectively ignored based on a logical mask.
%
% The three methods are:
% 1. Standard: Explicitly zeroing out rows/columns of A and using A*v.
% 2. Sparse Diagonal: The user-provided method using a sparse diagonal matrix.
% 3. Optimized Indexing: The proposed optimization using logical indexing.

clear;
clc;
close all;

% --- Configuration ---
N = 10000; % Matrix size (N x N). A larger N makes performance differences clearer.
fprintf('Starting benchmark with matrix size N = %d.\n', N);

% --- Initialization ---
% We will iterate from M=0 (no missing entries) to M=N (all missing).
% To avoid plotting all N+1 points which can be slow, we can take steps.
step_size = max(1, floor(N / 250)); % Take about 100 steps
missing_counts = [5000:step_size:8500];
num_steps = length(missing_counts);

% Pre-allocate arrays to store the timing results from timeit.
times_standard = zeros(1, num_steps);
times_sparse = zeros(1, num_steps);
times_optimized = zeros(1, num_steps);

% --- Main Benchmark Loop ---
% We use 'parfor' for parallel execution to speed up the benchmark process.
% If you don't have the Parallel Computing Toolbox, change 'parfor' to 'for'.
for i = 1:num_steps
    M = missing_counts(i);
    fprintf('Benchmarking with M = %d missing entries...\n', M);

    % --- Test Data Setup for this iteration ---
    % Create full random matrix and vector for consistent testing.
    A_full = rand(N, N);
    v = rand(N, 1);

    % Create the logical mask.
    mask = true(N, 1);
    if M > 0
        missing_indices = randperm(N, M);
        mask(missing_indices) = false;
    end

    % --- Method 1: Standard Multiplication ---
    % Create the explicitly zeroed-out matrix.
    A_masked = A_full;
    A_masked(~mask, :) = 0;
    A_masked(:, ~mask) = 0;
    % Time the standard multiplication.
    f1 = @() A_masked * v;
    times_standard(i) = timeit(f1);

    % --- Method 2: User's Sparse Diagonal Method ---
    f2 = @() masked_matvec_sparse(A_full, v, mask);
    times_sparse(i) = timeit(f2);

    % --- Method 3: Optimized Logical Indexing Method ---
    f3 = @() masked_matvec_optimized(A_full, v, mask);
    times_optimized(i) = timeit(f3);
end

fprintf('Benchmark finished.\n');

% --- Plotting Results ---
fprintf('Plotting results...\n');
figure('Name', 'Benchmark Results', 'NumberTitle', 'off', 'Color', 'w');
hold on;
grid on;

plot(missing_counts, times_standard * 1e3, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(missing_counts, times_sparse * 1e3, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(missing_counts, times_optimized * 1e3, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0, 0.6, 0.2]);

% --- Formatting the Plot ---
title(['Performance of Masked Matrix-Vector Multiplication (N = ' num2str(N) ')'], 'FontSize', 14);
xlabel('Number of Masked (Missing) Entries (M)', 'FontSize', 12);
ylabel('Execution Time (milliseconds)', 'FontSize', 12);
legend({'1. Standard (A_{masked} * v)', '2. Sparse Diagonal (Your Method)', '3. Optimized Indexing (My Method)'}, 'Location', 'northwest', 'FontSize', 10);
set(gca, 'YScale', 'log'); % Use a log scale for Y-axis to see all results clearly.
set(gca, 'FontSize', 10);
hold off;

% --- Calculate and Display Mean Times ---
mean_time_standard = mean(times_standard) * 1e3; % convert to ms
mean_time_sparse = mean(times_sparse) * 1e3;     % convert to ms
mean_time_optimized = mean(times_optimized) * 1e3; % convert to ms

fprintf('\n--- Arithmetic Mean of Execution Times ---\n');
fprintf('1. Standard Method:      %.4f ms\n', mean_time_standard);
fprintf('2. Sparse Diagonal Method: %.4f ms\n', mean_time_sparse);
fprintf('3. Optimized Indexing Method: %.4f ms\n', mean_time_optimized);
fprintf('------------------------------------------\n');


fprintf('Done.\n');


% --- Function Definitions ---
% These are included here to make the script self-contained.

function y = masked_matvec_sparse(A, x, mask)
    % User's original method.
    % A    : n×n double
    % mask : n×1 logical
    % x    : n×1 vector
    n = size(A,1);
    D = spdiags(mask, 0, n, n);
    y = D * (A * (D * x));
end

function y = masked_matvec_optimized(A, x, mask)
    % Optimized method using logical indexing.
    % A    : n×n double
    % mask : n×1 logical
    % x    : n×1 vector
    y = zeros(size(x), 'like', x);
    y(mask) = A(mask, mask) * x(mask);
end
