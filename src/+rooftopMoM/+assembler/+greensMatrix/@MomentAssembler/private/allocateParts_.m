function parts = allocateParts_(rooftops)
% allocateParts_: Returns the zero generators of the grid of the rooftops as a struct with the fields electric,
%   magneticXX and magneticYY; the generator of a missing rooftop direction is empty.

    primitive = rooftops.primitive();
    counts = [primitive.countX(), primitive.countY()];
    parts = struct('electric', zeros(counts, 'like', 1i), 'magneticXX', [], 'magneticYY', []);

    if rooftops.countDir('x') > 0
        parts.magneticXX = zeros(counts - [1, 0], 'like', 1i);
    end

    if rooftops.countDir('y') > 0
        parts.magneticYY = zeros(counts - [0, 1], 'like', 1i);
    end
end
