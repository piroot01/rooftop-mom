function circulant = generateTChan(generator, even)
% generateTChan: Returns the optimal circulant approximation of T. Chan of a Toeplitz generator, the circulant closest
%   in the Frobenius norm: every offset p, weighted by how often it occurs in the matrix, (count - |p|) / count, and
%   wrapped onto the period.
%
%   Inputs:
%       generator: (2 * countX - 1) x (2 * countY - 1) two-sided generator, or with even the countX x countY
%           one-sided generator of the offsets 0, ..., count - 1 (may be empty).
%       even: Optional true if the generator is one-sided and even, g(-p, q) = g(p, -q) = g(p, q) (default false).
%
%   Output:
%       circulant: countX x countY generator of the circulant.

    if isempty(generator)
        circulant = [];
        return;
    end

    circulant = rooftopMoM.assembler.preconditioner.utils.foldGenerator(generator, @window_, nargin > 1 && even);
end

function weights = window_(offsets, count, ~)
    weights = (count - abs(offsets)) / count;
end
