function spectralMatrixData = generateSplitSpectralGreensMatrix(magneticXX, magneticYY, electric, rooftopContent, faceSizeX, faceSizeY)
% generateSplitSpectralGreensMatrix: Generates the subband spectra of the Green's matrix for the split-FFT
%   accelerated product. A Toeplitz block T with a generator on n points is embedded into a circulant with the
%   period 2n, whose even and odd frequency subbands P and A give
%       T * x = (ifft(P .* fft(x)) + conj(w) .* ifft(A .* fft(w .* x))) / 2,  w(d) = exp(-1i * pi * d / n),
%   so the product needs only FFTs of the grid size. The splitting is done in every direction with more than one
%   face and n is rounded up to a fast FFT size.
%
%   Inputs:
%       magneticXX: Toeplitz generator of the magnetic part of the x-directed rooftops (empty for YONLY).
%       magneticYY: Toeplitz generator of the magnetic part of the y-directed rooftops (empty for XONLY).
%       electric: Toeplitz generator of the electric part spanning the whole face grid.
%       rooftopContent: RooftopContent of the rooftop basis.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Output:
%       spectralMatrixData: 1xS array of SpectralMatrixData with the subband blocks including the factor 1/2 per
%           split direction. For BOTH S = 4 and the element k = 1 + parityX + 2 * parityY holds the subband with
%           the given parities, otherwise S = 2 and the elements hold the even and the odd subband.

    [countX, countY] = size(electric);

    fftSizeX = 1;
    if countX > 1
        fftSizeX = rooftopMoM.core.greensMatrix.utils.generateFastFFTSize(countX);
    end

    fftSizeY = 1;
    if countY > 1
        fftSizeY = rooftopMoM.core.greensMatrix.utils.generateFastFFTSize(countY);
    end

    [paritiesX, paritiesY] = ndgrid(0:double(fftSizeX > 1), 0:double(fftSizeY > 1));
    scale = 1 / numel(paritiesX);

    for k = numel(paritiesX):-1:1
        [forwardShiftX, ~, forwardShiftY, ~] = rooftopMoM.core.greensMatrix.utils.generateSpectralShiftVectors( ...
            fftSizeX, fftSizeY, paritiesX(k), paritiesY(k));

        spectralMatrixData(k) = rooftopMoM.core.greensMatrix.utils.generateSpectralMatrixData( ...
            scale * rooftopMoM.core.greensMatrix.utils.generateSubbandSpectrum(magneticXX, fftSizeX, fftSizeY, paritiesX(k), paritiesY(k)), ...
            scale * rooftopMoM.core.greensMatrix.utils.generateSubbandSpectrum(magneticYY, fftSizeX, fftSizeY, paritiesX(k), paritiesY(k)), ...
            scale * rooftopMoM.core.greensMatrix.utils.generateSubbandSpectrum(electric, fftSizeX, fftSizeY, paritiesX(k), paritiesY(k)), ...
            rooftopContent, faceSizeX, faceSizeY, forwardShiftX, forwardShiftY ...
        );
    end
end
