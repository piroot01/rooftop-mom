function circulant = generateRChan(generator, even)
% generateRChan: Returns the R. Chan circulant approximation of a Toeplitz generator: every offset wrapped onto the
%   period, so the offset k > 0 gets the generator at k plus the generator at k - count.
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

    circulant = rooftopMoM.assembler.preconditioner.utils.foldGenerator(generator, @(offsets, ~, ~) ...
        ones(size(offsets)), nargin > 1 && even);
end
