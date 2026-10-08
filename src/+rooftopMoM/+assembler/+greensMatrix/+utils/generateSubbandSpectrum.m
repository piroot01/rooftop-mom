function spectrum = generateSubbandSpectrum(kernel, fftSizeX, fftSizeY, parityX, parityY)
% generateSubbandSpectrum: Returns one frequency subband of the spectrum of the circulant embedding (with the periods
%   2 * fftSizeX and 2 * fftSizeY) of a Toeplitz generator, without forming the embedding.
%
%   Inputs:
%       kernel: 2D Toeplitz generator (may be empty).
%       fftSizeX: Positive integer specifying the subband size in the x direction.
%       fftSizeY: Positive integer specifying the subband size in the y direction.
%       parityX: 0 for the even, 1 for the odd subband in the x direction (0 if fftSizeX is 1).
%       parityY: 0 for the even, 1 for the odd subband in the y direction (0 if fftSizeY is 1).
%
%   Output:
%       spectrum: fftSizeX x fftSizeY subband spectrum, empty for an empty kernel and real for a real kernel.

    if isempty(kernel)
        spectrum = [];
        return;
    end

    subbandKernel = rooftopMoM.assembler.greensMatrix.utils.generateSubbandKernel(kernel, fftSizeX, parityX, 1);
    subbandKernel = rooftopMoM.assembler.greensMatrix.utils.generateSubbandKernel(subbandKernel, fftSizeY, parityY, 2);

    spectrum = fft2(subbandKernel);

    % The spectrum of a real even circulant is real; drop the rounding in the imaginary part.
    if isreal(kernel)
        spectrum = real(spectrum);
    end
end
