function data = splitFFTMultDataPrep(toeplitzKernel)
% splitFFTMultDataPrep: Prepare necessary parts for the splitFFTMult algorithm.
%
%   Input:
%       toeplitzKernel: M-by-N Toeplitz kernel (numeric or complex). M, N must be >= 1.
%
%   Output (struct) data contains:
%       M, N: Size of the Toeplitz kernel.
%       hc_ee, hc_eo, hc_oe, hc_oo: Subsampled 2D FFT blocks of the circulant kernel.
%       pr, pc: Phase shifters.

    [M, N] = size(toeplitzKernel);
    data.M = M; data.N = N;

    % --- build circulant embedding hc (2M x 2N), same layout as original ---
    % Use 'like' to preserve class (single/double/complex).
    circulantKernel = zeros(2*M, 2*N, 'like', toeplitzKernel);

    % Top-left: original kernel.
    circulantKernel(1:M, 1:N) = toeplitzKernel;

    % Bottom-left: vertically flipped rows 2:end of toeplitzKernel.
    if M > 1
        circulantKernel((M + 2):end, 1:N) = flipud(toeplitzKernel(2:end, :));
    end

    % Top-right: horizontally flipped columns 2:end of toeplitzKernel.
    if N > 1
        circulantKernel(1:M, (N + 2):end) = fliplr(toeplitzKernel(:, 2:end));
    end

    % Bottom-right: rotate the interior 180° (rows 2:end, cols 2:end).
    if M > 1 && N > 1
        circulantKernel((M + 2):end, (N + 2):end) = rot90(toeplitzKernel(2:end, 2:end), 2);
    end

    % NOTE: Row M+1 and column N+1 remain zeros.

    % Compute 2D FFT and subsample the interleaved spectra.
    circulantKernelFFT = fft2(circulantKernel);

    data.circulantKernelFFT_ee = circulantKernelFFT(1:2:end, 1:2:end);
    data.circulantKernelFFT_eo = circulantKernelFFT(1:2:end, 2:2:end);
    data.circulantKernelFFT_oe = circulantKernelFFT(2:2:end, 1:2:end);
    data.circulantKernelFFT_oo = circulantKernelFFT(2:2:end, 2:2:end);

    % Build phase shift matrices used in split FFT multiply.
    krow = (0:(M-1)).';
    kcol = (0:(N-1)).';

    phaseShiftRow = exp(1i * pi * krow ./ M);
    phaseShiftColumn = exp(1i * pi * kcol ./ N);

    % Expand to full M-by-N matrices (explicit broadcasting).
    data.phaseShiftRow = repmat(phaseShiftRow, 1, N);
    data.phaseShiftColumn = repmat(phaseShiftColumn.', M, 1);
end
