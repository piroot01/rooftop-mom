clc; clearvars; close all;

repeatCount = 3;
Ms = round(logspace(log10(2), log10(1000), 100));
sizes = Ms.^2;
MsSize = numel(Ms);

timesOriginal = zeros(1, MsSize);
timesSplit = zeros(1, MsSize);
speedRatio = zeros(1, MsSize);

for n = 1:numel(Ms)
    M = Ms(n);
    fprintf('Computing with a vector size of %d, ', sizes(n));

    Ht = generateToeplitzKernel(M, M);
    X = generateRandomMatrix(M, M);
    x = X(:);

    localTimesOriginal = zeros(1, repeatCount);
    localTimesSplit = zeros(1, repeatCount);

    for r = 1:repeatCount
        [localTimesOriginal(r), localTimesSplit(n), normDiff] = runSingleComparison(Ht, x);
    end

    timesOriginal(n) = sum(localTimesOriginal) / repeatCount;
    timesSplit(n) = sum(localTimesSplit) / repeatCount;
    speedRatio(n) = timesOriginal(n) / timesSplit(n);
    fprintf('normDiff = %g.\n', normDiff);
end

% display the data
figure;
loglog(sizes, timesOriginal);
hold on;
loglog(sizes, timesSplit);

pTimesOriginal = polyfit(sizes, timesOriginal, 1);
pTimesSplit = polyfit(sizes, timesSplit, 1);

figure;
loglog(sizes, pTimesOriginal(1) * sizes + pTimesOriginal(2));
hold on;
loglog(sizes, pTimesSplit(1) * sizes + pTimesSplit(2));

figure;
semilogx(sizes, speedRatio);

% helper functions
function [timeOriginal, timeSplit, normDiff] = runSingleComparison(toeplitzKernel, x)
    % original FFT multiplication algorithm
    data = originalFFTMultDataPrep(toeplitzKernel);
    
    t0 = tic;
    bOriginal = originalFFTMult(data, x);
    timeOriginal = toc(t0);

    % split FFT multiplication algorithm
    data = splitFFTMultDataPrep(toeplitzKernel);
    
    t1 = tic;
    bSplit = splitFFTMult(data, x);
    timeSplit = toc(t1);

    normDiff = norm(bOriginal - bSplit);
end

function Ht = generateToeplitzKernel(M, N)
    Ht = generateRandomMatrix(M, N);
end

function X = generateRandomMatrix(M, N)
    X = randi(100, [M, N]) + 1i * randi(100, [M, N]);
end