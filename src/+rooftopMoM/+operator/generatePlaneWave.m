function excitation = generatePlaneWave(rooftops, wavenumber, direction, polarization, intensity)
% generatePlaneWave: Returns the excitation vector of an incident plane wave,
%   E(r) = intensity * polarization * exp(-1i * wavenumber * direction.' * r), tested with every rooftop.
%
%   Inputs:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%       wavenumber: Wavenumber of the wave.
%       direction: 3x1 unit vector, the direction of propagation.
%       polarization: 3x1 unit vector of the electric field, perpendicular to the direction.
%       intensity: Amplitude of the electric field in V/m.
%
%   Output:
%       excitation: Nx1 vector, the x-directed rooftops first.
%
%   NOTE: The rooftops of one direction are translates of the reference rooftop. Each entry is the entry of the
%       reference rooftop times the phase of the wave at the translation, so the cost is O(N).

    primitive = rooftops.primitive();
    faceArea = primitive.sizeX() * primitive.sizeY();
    [pointsU, pointsV] = rooftops.quadrature().points();
    weights = rooftops.quadrature().weights2D();

    % Exponent of the phase of the wave per face along x and y.
    phaseSteps = -1i * wavenumber * direction(1:2) .* [primitive.sizeX(); primitive.sizeY()];
    pointPhases = exp(phaseSteps(1) * pointsU + phaseSteps(2) * pointsV);

    directions = 'xy';
    parts = cell(2, 1);

    for n = 1:2
        rooftop = rooftops.referenceRooftop(directions(n));
        [rowCount, columnCount] = rooftops.gridSize(directions(n));

        secondFacePhase = exp(phaseSteps.' * rooftop.secondFaceOffset);
        amplitudes = rooftop.amplitudes(:, 1) + secondFacePhase * rooftop.amplitudes(:, 2);
        reference = intensity * polarization(n) * faceArea * pointPhases * (weights .* amplitudes);

        translations = exp(phaseSteps(1) * (0:(rowCount - 1)).') * exp(phaseSteps(2) * (0:(columnCount - 1)));
        parts{n} = reference * translations(:);
    end

    excitation = cat(1, parts{:});
end
