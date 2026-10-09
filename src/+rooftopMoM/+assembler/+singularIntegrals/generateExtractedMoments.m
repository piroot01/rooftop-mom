function moments = generateExtractedMoments(sizeX, sizeY, wavenumber)
% generateExtractedMoments: Returns the seven moments of the extracted part of the Green's function,
%   S(R) = 1 / R - k^2 R / 2 + k^4 R^3 / 24, over the cell pairs at the touching offsets (0, 0), (1, 0), (0, 1) and
%   (1, 1), in closed form. S holds the first three terms of the Taylor series of exp(-ikR) / R that are not smooth
%   at R = 0. On the touching cell pairs the quadrature integrates the remainder exp(-ikR) / R - S (evaluateRemainder),
%   whose derivatives up to the fourth order are continuous (the near-field note in docs/nearfield, Extraction of
%   order N and Regularity of the remainder).
%
%   Inputs:
%       sizeX, sizeY: Cell sizes.
%       wavenumber: Wavenumber k.
%
%   Output:
%       moments: 7x4 moments, ordered as by generatePowerMoments.

    moments = rooftopMoM.assembler.singularIntegrals.generatePowerMoments(sizeX, sizeY, -1) ...
        - wavenumber^2 / 2 * rooftopMoM.assembler.singularIntegrals.generatePowerMoments(sizeX, sizeY, 1) ...
        + wavenumber^4 / 24 * rooftopMoM.assembler.singularIntegrals.generatePowerMoments(sizeX, sizeY, 3);
end
