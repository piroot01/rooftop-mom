function result = splitFFTMult(data, input)
% splitFFTMult: Fast 2D multiply via split-FFT tree.
%
%   Inputs:
%       data: Struct produced by splitFFTMultDataPrep with fields:
%           M, N: Integer sizes (rows, cols).
%           phaseShiftRow: Row phase-shift tile (M-by-N).
%           phaseShiftColumn: Column phase-shift tile (M-by-N).
%           circulantKernelFFT_ee: Subsampled frequency-domain kernel (and _eo, _oe, _oo).
%       input: Column vector of length M*N (will be reshaped to M-by-N).
%
%   Output:
%       result: Column vector of length M*N (result of the split-FFT multiply).

    M = data.M;
    N = data.N;

    % Pull out the phase shift multipliers.
    phaseShiftRow = data.phaseShiftRow;
    phaseShiftColumn = data.phaseShiftColumn;

    % Reshape input vector into matrix (M rows, N cols).
    inputReshaped = reshape(input, M, N);

    % Optimise the tree structure based on the kernel size.
    if M >= N
        d1 = 1; d2 = 2; p1 = phaseShiftRow; p2 = phaseShiftColumn;
    else
        d1 = 2; d2 = 1; p1 = phaseShiftColumn; p2 = phaseShiftRow;
    end

    % Helper oneliners.
    % fft along dimension dim (preserve shape)
    fftd = @(input, dim) fft(input, [], dim);
    % Multiply by shift (broadcasting compatible) then fft and circshift -1 along dim.
    fftdShift = @(input, shift, dim) circshift(fft(input .* shift, [], dim), -1, dim);

    % Build FFT "tree" (expand)
    % First level (along d1).
    inputFFT_e  = fftd(inputReshaped, d1);
    inputFFT_o  = fftdShift(inputReshaped, p1, d1);

    % Second level (along d2): expand each of the two branches into two.
    inputFFT_ee = fftd(inputFFT_e, d2);
    inputFFT_eo = fftdShift(inputFFT_e, p2, d2);
    inputFFT_oe = fftd(inputFFT_o, d2);
    inputFFT_oo = fftdShift(inputFFT_o, p2, d2);

    % When M < N the eo and oe are swapped at this point.
    if M < N
        tmp = inputFFT_eo;
        inputFFT_eo = inputFFT_oe;
        inputFFT_oe = tmp;
    end

    % Pointwise multiplication in spectra.
    outputFFT_ee = data.circulantKernelFFT_ee .* inputFFT_ee;
    outputFFT_eo = data.circulantKernelFFT_eo .* inputFFT_eo;
    outputFFT_oe = data.circulantKernelFFT_oe .* inputFFT_oe;
    outputFFT_oo = data.circulantKernelFFT_oo .* inputFFT_oo;

    % When M < N the eo and oe are swapped at this point.
    if M < N
        tmp = outputFFT_eo;
        outputFFT_eo = outputFFT_oe;
        outputFFT_oe = tmp;
    end

    % Collapse the FFT tree (inverse transforms)
    % Second level inverse (along d2).
    outputFFT_e = ifft(outputFFT_ee, [], d2) + ifft(outputFFT_eo, [], d2) .* p2;
    outputFFT_o = ifft(outputFFT_oe, [], d2) + ifft(outputFFT_oo, [], d2) .* p2;

    % First level inverse (along d1).
    output = 0.25 .* (ifft(outputFFT_e, [], d1) + ifft(outputFFT_o, [], d1) .* p1);

    % Return as column vector.
    result = output(:);
end
