%% toeplitz_split_recursive.m
close all; clearvars;
% General d-dimensional Toeplitz matvec using recursive split-FFT (vectorized).
% Contains:
% - pre = toeplitz_precompute(ht)
% - b = toeplitz_matvec(pre, x)
% - demo at the bottom to compare with classical embedding

function plot_toeplitz_timings(startElems, endElems, numPoints, repeatCount)
% plot_toeplitz_timings  Run timing experiment and plot results.
%   plot_toeplitz_timings(startElems, endElems, numPoints, repeatCount)
%
% Inputs:
%   startElems  - starting number of elements (e.g. 1e4)
%   endElems    - ending number of elements (e.g. 9e6)
%   numPoints   - number of sizes to test (e.g. 8)
%   repeatCount - repeat each matvec this many times to average (e.g. 5)
%
% Notes:
%  - The function tries to pick M and N that are as close to square as possible
%    such that M*N >= targetElements (but not necessarily perfect square).
%  - If toeplitz_precompute / toeplitz_matvec are not found, the script
%    falls back to the classical FFT method (so the split-method timings
%    will match the classical timings unless you provide your implementations).

if nargin < 4, repeatCount = 5; end
if nargin < 3, numPoints = 6; end
if nargin < 2, endElems = 3e6; end
if nargin < 1, startElems = 1e4; end

% Use geometric spacing (logspace) so we cover wide ranges well
targets = round(logspace(log10(startElems), log10(endElems), numPoints));
% Make unique and sorted
targets = unique(max(2, targets));

rng(0); % reproducible random numbers

numSizes = numel(targets);
elems = zeros(numSizes,1);
time0 = zeros(numSizes,1);      % classical per-matvec time
time1 = zeros(numSizes,1);      % split-method per-matvec time
time0Prep = zeros(numSizes,1);  % classical prep time
time1Prep = zeros(numSizes,1);  % split-method prep time

fprintf('Running %d tests (repeatCount = %d)\n', numSizes, repeatCount);

for ii = 1:numSizes
    target = targets(ii);
    % Find M,N as close to square as possible where M*N >= target
    M = floor(sqrt(target));
    if M < 2, M = 2; end
    N = ceil(target / M);
    % try adjust M/N to reduce difference M*N - target but keep closeness to square
    % small local search (few steps)
    bestDiff = M*N - target;
    bestM = M; bestN = N;
    for d = -5:5
        Mcand = max(2, M + d);
        Ncand = ceil(target / Mcand);
        diff = Mcand * Ncand - target;
        % prefer smaller diff, then prefer Mcand closer to sqrt
        if diff >= 0
            if diff < bestDiff || (diff == bestDiff && abs(Mcand - Ncand) < abs(bestM - bestN))
                bestDiff = diff;
                bestM = Mcand; bestN = Ncand;
            end
        end
    end
    M = bestM; N = bestN;
    elems(ii) = M * N;
    fprintf('Test %d/%d : target=%d -> M=%d, N=%d (elements=%d)\n', ii, numSizes, target, M, N, elems(ii));
    
    % Prepare data
    ht = randi(100, [M, N]);
    xvec = randi(100, [M*N, 1]);
    
    % ----- Classical prep (mirrored circulant embedding) -----
    t0PrepStart = tic;
    rIdx = [1:M, (M - 1):-1:2];
    cIdx = [1:N, (N - 1):-1:2];
    hc0 = ht(rIdx, cIdx);
    hc_fft = fft2(hc0);
    time0Prep(ii) = toc(t0PrepStart);
    
    % ----- Classical matvec timing (method 1) -----
    t0start = tic;
    % repeat to average
    for k = 1:repeatCount
        x_reshaped = reshape(xvec, M, N);
        x_fft = fft2(x_reshaped, 2 * M - 2, 2 * N - 2);
        b_raw = ifft2(hc_fft .* x_fft);
        b_cutted = b_raw(1:M, 1:N);
        b0 = b_cutted(:);
    end
    time0(ii) = toc(t0start) / repeatCount;
    
    % ----- Split/precompute method prep -----
    t1PrepStart = tic;
    pre = split_fft_mult_prep_2d(ht);
    time1Prep(ii) = toc(t1PrepStart);
    
    % ----- Split matvec timing -----
    t1start = tic;
    for k = 1:repeatCount
        %b1 = split_fft_mult(pre, xvec);
        b1 = split_fft_mult_2d(pre, xvec);
    end
    time1(ii) = toc(t1start) / repeatCount;
    
    % Small pause to let system breathe (optional, comment out if you want speed)
    pause(0.01);
end

% -------------------- Plotting --------------------
figure('Color','w','Position',[200 200 900 560]);
% choose line widths and markers for readability
lw = 1.8;
ms = 8;

% Plot times
hold on;
% classical & split matvec times (solid)
h1 = plot(elems, time0, '-o', 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Classical matvec (time0)');
h2 = plot(elems, time1, '-s', 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Split matvec (time1)');

% precomputation times (dashed)
h3 = plot(elems, time0Prep, '--o', 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Classical prep (time0Prep)');
h4 = plot(elems, time1Prep, '--s', 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Split prep (time1Prep)');

%plot(time0 / time1, 'LineWidth', 2);

% Aesthetics
grid on;
hold off;

% Print a summary table
T = table((1:numSizes)', elems, time0, time0Prep, time1, time1Prep, ...
    'VariableNames', {'idx','elements','time0_s','time0Prep_s','time1_s','time1Prep_s'});
disp('Summary results:');
disp(T);

end

plot_toeplitz_timings(1e4, 9e6, 10, 10);
