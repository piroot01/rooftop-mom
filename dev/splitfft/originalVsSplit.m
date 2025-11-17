%% FFT vs split FFT matrix-vector multiplication benchmark
%% Single-run comparison
% This script compares two approaches for accelerating a matrix-vector
% multiplication that arises from a 2D Toeplitz-like convolution:
%  - FFTMult: a straightforward 2-D FFT-based approach
%  - splitFFTMult: a memory- and performance-optimized approach that
%    splits the 2-D FFT into row/column operations and applies phase
%    shifts to reuse computations.
clc; clearvars; close all;

M = 100;
N = M;

ht = complex(randn(M, N), randn(M, N));
x = complex(randn(M, N), randn(M, N));

% --------------------------------------------------
% 2-D FFT-accelerated matrix-vector multiplication
% --------------------------------------------------
function data = FFTMultDataPrep_(ht)
    [M, N] = size(ht);

    % allocate zero-padded, mirrored kernel
    hc = zeros(2*M, 2*N, 'like', ht);
    hc(1:M, 1:N) = ht;                                  % top-left block
    hc((M + 2):end, 1:N) = flipud(ht(2:end, :));       % bottom-left mirrored rows
    hc(1:M, (N + 2):end) = fliplr(ht(:, 2:end));       % top-right mirrored cols
    hc((M + 2):end, (N + 2):end) = rot90(ht(2:end, 2:end), 2); % bottom-right rotated block

    % store the Fourier transform and sizes for later use
    data = struct( ...
        'hcf', fft2(hc), ...
        'M', M, ...
        'N', N ...
    );
end

function y = FFTMult_(data, x)
    M = data.M;
    N = data.N;

    % compute 2-D FFT of the input (zero-padded to 2M-by-2N) and
    % perform pointwise multiplication in the frequency domain
    xf = fft2(x, 2*M, 2*N);
    yf = data.hcf .* xf;

    % free the intermediate to reduce peak memory usage
    clear xf;

    % inverse FFT and crop back to the original domain
    y = ifft2(yf);
    y = y(1:M, 1:N);
end

function memoryB = predictMemoryUsageFFTMult(M, N, memoryOverheadB, elementMemoryB)
    % Rough estimate of peak memory usage (in bytes) for the
    % straightforward 2-D FFT approach.

    if nargin < 4
        elementMemoryB = 16; % assume complex double: 16 bytes per element
    end

    if nargin < 3
        % default overhead (safety margin) set to 1 GB
        memoryOverheadB = 1024^3;
    end

    % memory used by ht (original kernel)
    htMemoryB = elementMemoryB * M * N;

    % memory used by x (input) - same size as ht
    xMemoryB = htMemoryB;

    % hc (zero-padded mirror) is 4 times larger than ht
    hcMemoryB = 4 * htMemoryB;
    
    % account for additional FFT allocations: two more copies of hc
    % (this is conservative) and include overhead
    memoryB = htMemoryB + xMemoryB + 3 * hcMemoryB + memoryOverheadB;
end

% --------------------------------------------------
% 2-D split-FFT accelerated matrix-vector multiplication
% --------------------------------------------------
function data = splitFFTMultDataPrep_(ht)
    data = FFTMultDataPrep_(ht);

    M = data.M;
    N = data.N;

    % extract sub-sampled frequency-domain blocks: these correspond to
    % the 2x2 splitting of the padded transform (even/odd rows and cols)
    data.hceef = data.hcf(1:2:end, 1:2:end); % even rows, even cols
    data.hceof = data.hcf(1:2:end, 2:2:end); % even rows, odd cols
    data.hcoef = data.hcf(2:2:end, 1:2:end); % odd rows, even cols
    data.hcoof = data.hcf(2:2:end, 2:2:end); % odd rows, odd cols

    % free the full transform to save memory
    data.hcf = [];
    
    % build row/column index vectors for phase multipliers
    kr = (0:(M - 1)).';
    kc = (0:(N - 1));

    % phase factors used to shift between even/odd split domains
    data.pr = exp(-1i * pi * kr / M);
    data.pc = exp(-1i * pi * kc / N);
    data.cpr = conj(data.pr);
    data.cpc = conj(data.pc);
end

function y = splitFFTMult_(data, x)
    % compute FFT across rows (dimension 1)
    X_rows = fft(x, [], 1);

    % multiply with the appropriate frequency-domain submatrices and
    % perform column FFTs (dimension 2)
    Yee = data.hceef .* fft(X_rows, [], 2);     % even-even block
    Yeo = data.hceof .* fft(X_rows .* data.pc, [], 2);  % even-odd block

    % free row-FFT intermediates to reduce peak memory usage
    clear X_rows;

    % prepare row-FFT when input rows are multiplied by pr (phase shift)
    X_rows_pr = fft(x .* data.pr, [], 1);

    Yoe = data.hcoef .* fft(X_rows_pr, [], 2);              % odd-even block
    Yoo = data.hcoof .* fft(X_rows_pr .* data.pc, [], 2);   % odd-odd block

    clear X_rows_pr;

    % inverse FFT along columns and combine even/odd contributions
    % cpc contains conjugates of the column-phase factors
    Ye = ifft(Yee, [], 2) + ifft(Yeo, [], 2) .* data.cpc;
    clear Yee Yeo;

    Yo  = ifft(Yoe, [], 2) + ifft(Yoo, [], 2) .* data.cpc;
    clear Yoe Yoo;

    y = 0.25 * (ifft(Ye, [], 1) + ifft(Yo, [], 1) .* data.cpr);
end

function memoryB = predictMemoryUsageSplitFFTMult(M, N, memoryOverheadB, elementMemoryB)
    % Estimate peak memory usage (in bytes) for the split FFT method.

    if nargin < 4
        elementMemoryB = 16;
    end

    if nargin < 3
        % default overhead (safety margin) set to 1 GB
        memoryOverheadB = 1024^3;
    end

    % ht memory
    htMemoryB = elementMemoryB * M * N;

    % x memory
    xMemoryB = htMemoryB;

    % full padded hc memory (before subsampling) is 4 * ht
    hcMemoryB = 4 * htMemoryB;

    % memory for the phase shift vectors (pr, pc and their conjugates)
    pMemoryB = 2 * (M + N) * elementMemoryB;
    
    % combine the contributions. The formula attempts to capture the
    % memory used by input, kernel, phase vectors, and temporary
    % buffers used during FFTs. The constants are conservative.
    memoryB = (xMemoryB + htMemoryB + pMemoryB + hcMemoryB) + ((4 + 2) * xMemoryB) + memoryOverheadB;
end

% --------------------------------------------------
% Run quick sanity check and print memory estimates
% --------------------------------------------------
memoryB = predictMemoryUsageFFTMult(M, N, 0);
fprintf('FFTMult max memory requirements: %g GB\n', memoryB / 1024^3);

tic;
data = FFTMultDataPrep_(ht);
toc;

tic;
y_FFTMult = FFTMult_(data, x);
toc;

%clear data y_FFTMult;

memoryB = predictMemoryUsageSplitFFTMult(M, N, 0);
fprintf('splitFFTMult max memory requirements: %g GB\n', memoryB / 1024^3);

tic;
data = splitFFTMultDataPrep_(ht);
toc;

tic;
y_splitFFTMult = splitFFTMult_(data, x);
toc;

%clear data y_splitFFTMult;

norm(y_FFTMult - y_splitFFTMult)

%% memory comparison
clc; clearvars; close all;

BinGB = 1024^3;

% Get the free memory available on the system
%freeMemoryB = getMemoryInfo();
%freeMemoryGB = freeMemoryB / BinGB;

freeMemoryGB = 200;
freeMemoryB = freeMemoryGB * BinGB;

% inputs
count = 100;
Mmin = 100;
Mmax = 1e5;

Ms = round(logspace(log10(Mmin), log10(Mmax), count));
MsSq = Ms.^2;

% evaluate estimated memory usage for a range of square problems
memoryB_FFTMult = arrayfun(@(M) predictMemoryUsageFFTMult(M, M, 0), Ms);
memoryGB_FFTMult = memoryB_FFTMult ./ BinGB;
memoryB_splitFFTMult = arrayfun(@(M) predictMemoryUsageSplitFFTMult(M, M, 0), Ms);
memoryGB_splitFFTMult = memoryB_splitFFTMult ./ BinGB;

% find the largest problem that fits into available RAM
maxFFTMultIdx = find(memoryB_FFTMult < freeMemoryB, 1, 'last');
maxSplitFFTMultIdx = find(memoryB_splitFFTMult < freeMemoryB, 1, 'last');

fprintf('Free memory: %g GB.\n', freeMemoryGB);
fprintf('Maximum number of unknowns computable by FFTMult from the Ms: %g (M = %g).\n', MsSq(maxFFTMultIdx), Ms(maxFFTMultIdx));
fprintf('Maximum number of unknowns computable by splitFFTMult from the Ms: %g (M = %g).\n', MsSq(maxSplitFFTMultIdx), Ms(maxSplitFFTMultIdx));

figure('Name', 'memory comparison');
h1 = loglog(MsSq, memoryGB_FFTMult, 'LineWidth', 2);
hold on;
h2 = loglog(MsSq, memoryGB_splitFFTMult, 'LineWidth', 2);
hold on;
h3 = loglog(MsSq, ones(size(MsSq)) * freeMemoryGB, 'LineWidth', 2, 'LineStyle', '--');
hold off;
xlim([MsSq(1), MsSq(end)]);
xlabel('Unknown count [-]');
ylabel('Memory [GB]');
legend([h1 h2 h3], {'FFTMult','splitFFTMult','free memory'}, 'Location', 'best');

function [freeMemoryB, usedMemoryB, totalMemoryB] = getMemoryInfo()
    % Platform-dependent query for total and available physical memory.
    % Returns freeMemoryB (available), usedMemoryB, and totalMemoryB.
    if ispc
        [~, sys] = memory;
        freeMemoryB = sys.PhysicalMemory.Available;
    elseif ismac
        % macOS: use sysctl and vm_stat to estimate available memory
        [~, out1] = system('sysctl -n hw.memsize');
        totalMemoryB = str2double(strtrim(out1));

        % Free/available RAM: use vm_stat to get "Pages free" and "Pages inactive"
        [~, out2] = system('vm_stat');
        % page size is reported in the first line like "Mach Virtual Memory Statistics: (page size of 4096 bytes)"
        pg = regexp(out2, 'page size of (\d+)', 'tokens', 'once');
        pageSize = str2double(pg{1});
        % parse numeric page counts
        nums = regexp(out2, '([A-Za-z ]+):\s+(\d+)\.', 'tokens');
        freePages = 0;
        for k = 1:numel(nums)
            name = strtrim(nums{k}{1});
            count = str2double(nums{k}{2});
            if strcmpi(name, 'Pages free') || strcmpi(name, 'Pages inactive')
                freePages = freePages + count;
            end
        end
        freeMemoryB = freePages * pageSize;
    else
        % Linux: read /proc/meminfo for MemTotal and MemAvailable
        txt = fileread('/proc/meminfo');
        tokTotal = regexp(txt, '^MemTotal:\s+(\d+)\s+kB', 'tokens', 'lineanchors');
        tokAvail = regexp(txt, '^MemAvailable:\s+(\d+)\s+kB', 'tokens', 'lineanchors');
        totalMemoryB = str2double(tokTotal{1}) * 1024;
        freeMemoryB  = str2double(tokAvail{1})  * 1024;
    end

    usedMemoryB = totalMemoryB - freeMemoryB;
end

%% accuracy comparison
clc; clearvars; close all;

fprintf('\n[[ Running accuracy comparison ]]\n\n');

BinGB = 1024^3;

% Choose the memory
freeMemoryB = getMemoryInfo();
freeMemoryGB = freeMemoryB / BinGB;

% inputs
count = 10;
Mmin = 1e2;
Mmax = 1e3;

Ms = round(logspace(log10(Mmin), log10(Mmax), count));
MsSq = Ms.^2;

% check the memory requirements of this script
% we assume that the split FFT multiplication is less memory demanding and
% that the results fit in available RAM
maxMemoryB_FFTMult = predictMemoryUsageFFTMult(Mmax, Mmax);

fprintf('Maximal expected memory usage %g GB.\n', maxMemoryB_FFTMult / BinGB);

if freeMemoryB < maxMemoryB_FFTMult
    error('Not enough free memory.\n');
end

euclideanDistances = zeros(count, 1);

for n = 1:count
    fprintf('[%d/%d] Computing with %g (M = %g) unknowns.\n', n, count, MsSq(n), Ms(n));
    M = Ms(n);
    N = M;

    ht = complex(randn(M, N), randn(M, N));
    x = complex(randn(M, N), randn(M, N));

    y_FFTMult = FFTMult_(FFTMultDataPrep_(ht), x);
    y_splitFFTMult = splitFFTMult_(splitFFTMultDataPrep_(ht), x);
    euclideanDistances(n) = euclideanDistance(y_FFTMult(:), y_splitFFTMult(:));

    clear ht x y_FFTMult y_splitFFTMult;
end

figure('Name', 'accuracy comparison');
semilogx(MsSq, euclideanDistances, 'LineWidth', 2);

function distance = euclideanDistance(vec1, vec2)
    % Euclidean norm between two vectors
    distance = norm(vec1 - vec2);
end

%% speed comparison
clc; clearvars; close all;

fprintf('\n[[ Running speed comparison ]]\n\n');

BinGB = 1024^3;

% Choose the memory
freeMemoryB = getMemoryInfo();
freeMemoryGB = freeMemoryB / BinGB;

% inputs
warmupRunCount = 1;
runCount = 3;

% Choose the sample count, minimum and maximum problem sizes
count = 100;
Mmin = 1e2;
Mmax = 5e2;

enable_FFTMult = true;
enable_splitFFTMult = true;

Ms = round(logspace(log10(Mmin), log10(Mmax), count));
MsSq = Ms.^2;

timesS_FFTMult = nan(size(Ms));
timesS_splitFFTMult = nan(size(Ms));

for n = 1:count
    fprintf('[%d/%d] Computing with %g (M = %g) unknowns.\n', n, count, MsSq(n), Ms(n));
    M = Ms(n);
    N = M;
    
    if enable_FFTMult
        memoryB_FFTMult = predictMemoryUsageFFTMult(M, M, 0);
        if freeMemoryB < memoryB_FFTMult
            fprintf('Not enough free memory for the FFTMult (%g GB).\n', memoryB_FFTMult / BinGB);
            enable_FFTMult = false;
        end
    end

    if enable_splitFFTMult
        memoryB_splitFFTMult = predictMemoryUsageSplitFFTMult(M, M, 0);
        if freeMemoryB < memoryB_splitFFTMult
            fprintf('Not enough free memory for the splitFFTMult (%g GB).\n', memoryB_splitFFTMult / BinGB);
            enable_splitFFTMult = false;
        end
    end

    if ~enable_splitFFTMult & ~enable_splitFFTMult
        break;
    end

    ht = complex(randn(M, N), randn(M, N));
    x = complex(randn(M, N), randn(M, N));

    for m = 1:warmupRunCount
        if enable_FFTMult
            y_FFTMult = FFTMult_(FFTMultDataPrep_(ht), x);

            clear y_FFTMult;
        end

        if enable_splitFFTMult
            y_splitFFTMult = splitFFTMult_(splitFFTMultDataPrep_(ht), x);

            clear y_splitFFTMult;
        end
    end

    timeS_FFTMult = zeros(runCount, 1);
    timeS_splitFFTMult = zeros(runCount, 1);

    for m = 1:runCount
        if enable_FFTMult
            data = FFTMultDataPrep_(ht);
    
            t_FFTMult = tic;
            y_FFTMult = FFTMult_(FFTMultDataPrep_(ht), x);
            timeS_FFTMult(m) = toc(t_FFTMult);
    
            clear data y_FFTMult;
        end

        if enable_splitFFTMult
            data = splitFFTMultDataPrep_(ht);
    
            t_splitFFTMult = tic;
            y_splitFFTMult = splitFFTMult_(data, x);
            timeS_splitFFTMult(m) = toc(t_splitFFTMult);
    
            clear data y_splitFFTMult;
        end
        
        fprintf('\t[%d/%d] timeS_FFTMult = %g s, timeS_splitFFTMult = %g s, speed up %g\n', m, runCount, ...
            timeS_FFTMult(m), timeS_splitFFTMult(m), timeS_FFTMult(m) / timeS_splitFFTMult(m));
    end

    timesS_FFTMult(n) = mean(timeS_FFTMult);
    timesS_splitFFTMult(n) = mean(timeS_splitFFTMult);

    fprintf('[%d/%d] timesS_FFTMult = %g s, timesS_splitFFTMult = %g s, speed up %g\n', n, count, ...
            timesS_FFTMult(n), timesS_splitFFTMult(n), timesS_FFTMult(n) / timesS_splitFFTMult(n));

    clear ht x;
end

figure('Name', 'speed comparison');
h1 = loglog(MsSq, timesS_FFTMult, 'LineWidth', 2);
hold on;
h2 = loglog(MsSq, timesS_splitFFTMult, 'LineWidth', 2);
hold off;
xlim([MsSq(1), MsSq(end)]);
xlabel('Unknown count [-]');
ylabel('Time [s]');
legend([h1 h2], {'FFTMult','splitFFTMult'}, 'Location', 'best');
