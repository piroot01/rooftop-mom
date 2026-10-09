function circulant = generateRChan(generator)
% generateRChan: Returns the R. Chan circulant approximation of a Toeplitz generator: every offset wrapped onto the
%   period, so the offset k > 0 gets the generator at k plus the generator at k - count.
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

    circulant = rooftopMoM.assembler.preconditioner.utils.foldGenerator(generator, @(offsets, ~, ~) ...
        ones(size(offsets)));
end
