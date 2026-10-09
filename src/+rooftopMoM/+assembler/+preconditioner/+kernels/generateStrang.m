function circulant = generateStrang(generator)
% generateStrang: Returns the Strang circulant approximation of a Toeplitz generator: the offsets below half the
%   period, wrapped onto it. An offset of exactly half the period shares its place with its negative, so an even
%   generator gives a symmetric circulant.
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
    weights = double(abs(offsets) < count / 2) + 0.5 * (abs(offsets) == count / 2);
end
