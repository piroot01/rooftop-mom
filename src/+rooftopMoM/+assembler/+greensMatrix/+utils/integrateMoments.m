function moments = integrateMoments(rule, offsets, sizes, wavenumber)
% integrateMoments: Returns the seven moments of the Green's function exp(-ikR) / R over the cell pairs at the given
%   offsets, integrated with a rule of generateReducedRule. On the touching pairs, offsets in {0, 1}^2, the rule
%   integrates the remainder of the singularity extraction (singularIntegrals.evaluateRemainder) instead; the moments
%   of the extracted part are not included.
%
%   Inputs:
%       rule: Struct of generateReducedRule.
%       offsets: 2xC nonnegative cell offsets (dX; dY) of the other cell from the cell at the origin.
%       sizes: 1x2 cell sizes [sizeX, sizeY].
%       wavenumber: Wavenumber k.
%
%   Output:
%       moments: 7xC moments divided by the squared cell area, ordered as by singularIntegrals.generatePowerMoments.

    count = size(offsets, 2);

    % R = sqrt(a^2 (nu - dX)^2 + b^2 (eta - dY)^2), Px x Py x C. The singular points lie on the ends of the parts of
    % the rule, never on its nodes, so R > 0.
    distance = hypot(sizes(1) * (rule.nodesX - reshape(offsets(1, :), 1, 1, count)), ...
        sizes(2) * (rule.nodesY - reshape(offsets(2, :), 1, 1, count)));
    kernel = exp(-1i * wavenumber * distance) ./ distance;

    touching = find(all(offsets <= 1, 1));

    if ~isempty(touching)
        kernel(:, :, touching) = rooftopMoM.assembler.singularIntegrals.evaluateRemainder(distance(:, :, touching), ...
            wavenumber);
    end

    % Freed before the product, so the chunk peaks with the distances and the kernel, 24 bytes per kernel value.
    clear distance

    moments = rule.weights * reshape(kernel, [], count);
end
