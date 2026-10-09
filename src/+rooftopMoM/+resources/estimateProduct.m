function [stored, working] = estimateProduct(rooftops, product)
% estimateProduct: Returns the memory of a product with the Green's matrix before it is built: the bytes of its
%   spectra and the bytes alive at the peak of one product.
%
%   Inputs:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%       product: 'fft', 'splitFFT' or 'leanSplitFFT'.
%
%   Outputs:
%       stored: Bytes of the spectra.
%       working: Bytes of the arrays alive at the peak of one product.

    primitive = rooftops.primitive();
    counts = [primitive.countX(), primitive.countY()];
    content = rooftops.content();
    blockCount = 1 + 3 * (content == rooftopMoM.basisFunctions.RooftopContent.BOTH);

    if strcmp(product, 'fft')
        % The circulant embedding mirrors every direction with more than one face onto 2 * count - 2 points.
        sizes = max(2 * counts - 2, 1);
        subbandCount = 1;
    else
        % Every direction with more than one face splits into two subbands of a fast FFT size.
        sizes = [1, 1];

        for dimension = find(counts > 1)
            sizes(dimension) = rooftopMoM.assembler.greensMatrix.utils.generateFastFFTSize(counts(dimension));
        end

        subbandCount = 2^nnz(counts > 1);
    end

    applicators = struct('fft', 'FFTApplicator', 'splitFFT', 'SplitFFTApplicator', ...
        'leanSplitFFT', 'LeanSplitFFTApplicator');
    stored = 16 * blockCount * subbandCount * prod(sizes);
    working = feval(['rooftopMoM.assembler.greensMatrix.', applicators.(product), '.estimateWorkingMemory'], ...
        rooftops, sizes(1), sizes(2), content);
end
