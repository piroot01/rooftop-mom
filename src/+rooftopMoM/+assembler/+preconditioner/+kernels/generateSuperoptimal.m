function circulant = generateSuperoptimal(generator)
% generateSuperoptimal: Returns the superoptimal circulant approximation of Tyrtyshnikov of a Toeplitz generator, the
%   circulant C that minimizes norm(I - C \ T, 'fro'):
%       C = c(T * T') / c(T'),
%   with c the optimal circulant of T. Chan. c(T * T') follows from the parts of T without forming the product
%   (utils.generateOptimalOfProduct), so the cost stays that of a few FFTs.
%
%   Input:
%       generator: (2 * countX - 1) x (2 * countY - 1) two-sided generator (may be empty).
%
%   Output:
%       circulant: countX x countY generator of the circulant.

    if isempty(generator)
        circulant = [];
        return;
    end

    parts = rooftopMoM.assembler.preconditioner.utils.generatePartSpectra(generator);
    squared = rooftopMoM.assembler.preconditioner.utils.generateOptimalOfProduct(parts, parts);
    optimal = fft2(rooftopMoM.assembler.preconditioner.kernels.generateTChan(generator));

    % c(T') = c(T)' has the conjugate eigenvalues of c(T).
    circulant = ifft2(squared ./ conj(optimal));
end
