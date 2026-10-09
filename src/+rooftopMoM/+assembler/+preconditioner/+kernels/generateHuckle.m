function circulant = generateHuckle(generator, bandwidthX, bandwidthY)
% generateHuckle: Returns the Huckle circulant approximation of a Toeplitz generator: the offsets below the bandwidth,
%   damped linearly by 1 - |offset| / bandwidth in every direction and wrapped onto the period.
%
%   Inputs:
%       generator: (2 * countX - 1) x (2 * countY - 1) two-sided generator (may be empty).
%       bandwidthX: Positive integer specifying the bandwidth in the x direction.
%       bandwidthY: Positive integer specifying the bandwidth in the y direction.
%
%   Output:
%       circulant: countX x countY generator of the circulant.

    if isempty(generator)
        circulant = [];
        return;
    end

    bandwidths = [bandwidthX, bandwidthY];
    circulant = rooftopMoM.assembler.preconditioner.utils.foldGenerator(generator, ...
        @(offsets, ~, dimension) max(0, 1 - abs(offsets) / bandwidths(dimension)));
end
