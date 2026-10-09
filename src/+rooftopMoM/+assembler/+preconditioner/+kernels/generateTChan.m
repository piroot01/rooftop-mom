function circulant = generateTChan(generator)
% generateTChan: Returns the optimal circulant approximation of T. Chan of a Toeplitz generator, the circulant closest
%   in the Frobenius norm: every offset p, weighted by how often it occurs in the matrix, (count - |p|) / count, and
%   wrapped onto the period.
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

    circulant = rooftopMoM.assembler.preconditioner.utils.foldGenerator(generator, @window_);
end

function weights = window_(offsets, count, ~)
    weights = (count - abs(offsets)) / count;
end
