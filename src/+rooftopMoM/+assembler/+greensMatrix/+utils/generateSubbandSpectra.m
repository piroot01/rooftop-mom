function spectra = generateSubbandSpectra(generators, faceSizeX, faceSizeY)
% generateSubbandSpectra: Returns the subband block spectra of the Green's matrix for the split-FFT product.
%   A Toeplitz block T with a generator on n points is embedded into a circulant of the period 2n. Its even and odd
%   frequency subbands P and A give
%       T * x = (ifft(P .* fft(x)) + conj(w) .* ifft(A .* fft(w .* x))) / 2,  w(d) = exp(-1i * pi * d / n),
%   so the product needs FFTs of the grid size only. Every direction with more than one face is split, and n is
%   rounded up to a fast FFT size. One FFT of the embedding of a generator gives all its subbands: they are its even
%   and its odd frequencies along every split direction.
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

    start = tic;
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
    [stepX, stepY] = size(paritiesX);
    kernels = {generators.magneticXX(), generators.magneticYY(), generators.electric()};
    subbands = cell(3, numel(paritiesX));

    for n = find(~cellfun(@isempty, kernels))
        % The embedding of the period 2 * fftSize in every split direction holds t(p, q) at the offsets (+-p, +-q),
        % already divided by the number of subbands; the offset fftSize stays zero.
        [rows, columns] = size(kernels{n});
        embedding = zeros(stepX * fftSizeX, stepY * fftSizeY, 'like', kernels{n});
        embedding(1:rows, 1:columns) = kernels{n} / numel(paritiesX);
        embedding((end - rows + 2):end, 1:columns) = embedding(rows:(-1):2, 1:columns);
        embedding(:, (end - columns + 2):end) = embedding(:, columns:(-1):2);
        embedding = fft2(embedding);

        % The spectrum of a real even circulant is real; drop the rounding in the imaginary part.
        if isreal(kernels{n})
            embedding = real(embedding);
        end

        for k = 1:numel(paritiesX)
            subbands{n, k} = embedding((1 + paritiesX(k)):stepX:end, (1 + paritiesY(k)):stepY:end);
        end

        clear embedding
    end

    spectra = cell(1, numel(paritiesX));

    for k = 1:numel(paritiesX)
        [shiftX, shiftY] = rooftopMoM.assembler.greensMatrix.utils.generateShiftPhases(fftSizeX, fftSizeY, ...
            paritiesX(k), paritiesY(k));

        spectra{k} = rooftopMoM.assembler.greensMatrix.utils.generateBlockSpectra(subbands{:, k}, ...
            generators.content(), faceSizeX, faceSizeY, shiftX, shiftY);

        % Free the subbands of the generators as soon as their blocks exist, which keeps the peak memory low.
        subbands(:, k) = {[]};
    end

    spectra = [spectra{:}];

    rooftopMoM.logging.Logger.instance().info('Generated the subband spectra in %.3f s (%.1f MB).', toc(start), ...
        spectra.memory() / 2^20);
end
