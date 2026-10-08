function subbandKernel = generateSubbandKernel(kernel, fftSize, parity, dim)
% generateSubbandKernel: Folds a Toeplitz generator t into the kernel of one frequency subband of its circulant
%   embedding with the period 2 * fftSize. The even subband (parity 0) has the periodic kernel t(d) + t(fftSize - d),
%   the odd subband (parity 1) the modulated antiperiodic kernel (t(d) - t(fftSize - d)) .* exp(-1i * pi * d / fftSize).
%
%   Inputs:
%       kernel: 2D Toeplitz generator, its size in the direction dim must not exceed fftSize.
%       fftSize: Positive integer specifying the subband size in the direction dim.
%       parity: 0 for the even subband, 1 for the odd subband.
%       dim: Direction of the folding, 1 or 2.
%
%   Output:
%       subbandKernel: Subband kernel with the size fftSize in the direction dim.

    % Zero-pad the generator to fftSize + 1 samples, so t(fftSize) = 0.
    paddingSize = size(kernel);
    paddingSize(dim) = fftSize + 1 - paddingSize(dim);
    kernel = cat(dim, kernel, zeros(paddingSize, 'like', kernel));

    direct = {':', ':'};
    direct{dim} = 1:fftSize;

    reflected = {':', ':'};
    reflected{dim} = (fftSize + 1):(-1):2;

    if parity
        subbandKernel = (kernel(direct{:}) - kernel(reflected{:})) .* ...
            rooftopMoM.assembler.greensMatrix.utils.generateSubbandModulation(fftSize, dim);
    else
        subbandKernel = kernel(direct{:}) + kernel(reflected{:});
    end
end
