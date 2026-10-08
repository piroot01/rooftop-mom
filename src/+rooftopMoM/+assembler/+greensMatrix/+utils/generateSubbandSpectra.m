function spectra = generateSubbandSpectra(generators, faceSizeX, faceSizeY)
% generateSubbandSpectra: Returns the subband block spectra of the Green's matrix for the split-FFT product.
%   A Toeplitz block T with a generator on n points is embedded into a circulant of the period 2n. Its even and odd
%   frequency subbands P and A give
%       T * x = (ifft(P .* fft(x)) + conj(w) .* ifft(A .* fft(w .* x))) / 2,  w(d) = exp(-1i * pi * d / n),
%   so the product needs FFTs of the grid size only. Every direction with more than one face is split, and n is
%   rounded up to a fast FFT size.
%
%   Inputs:
%       generators: rooftopMoM.assembler.greensMatrix.Generators.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Output:
%       spectra: 1xS array of Spectra with the factor 1/2 per split direction included. With both directions split
%           S = 4 and the element 1 + parityX + 2 * parityY holds the subband of the given parities, otherwise S = 2
%           and the elements hold the even and the odd subband.

    [countX, countY] = size(generators.electric());

    fftSizeX = 1;
    if countX > 1
        fftSizeX = rooftopMoM.assembler.greensMatrix.utils.generateFastFFTSize(countX);
    end

    fftSizeY = 1;
    if countY > 1
        fftSizeY = rooftopMoM.assembler.greensMatrix.utils.generateFastFFTSize(countY);
    end

    [paritiesX, paritiesY] = ndgrid(0:double(fftSizeX > 1), 0:double(fftSizeY > 1));
    scale = 1 / numel(paritiesX);
    subband = @(kernel, k) scale * rooftopMoM.assembler.greensMatrix.utils.generateSubbandSpectrum(kernel, ...
        fftSizeX, fftSizeY, paritiesX(k), paritiesY(k));

    spectra = cell(1, numel(paritiesX));

    for k = 1:numel(paritiesX)
        [shiftX, shiftY] = rooftopMoM.assembler.greensMatrix.utils.generateShiftPhases(fftSizeX, fftSizeY, ...
            paritiesX(k), paritiesY(k));

        spectra{k} = rooftopMoM.assembler.greensMatrix.utils.generateBlockSpectra( ...
            subband(generators.magneticXX(), k), subband(generators.magneticYY(), k), ...
            subband(generators.electric(), k), generators.content(), faceSizeX, faceSizeY, shiftX, shiftY);
    end

    spectra = [spectra{:}];
end
