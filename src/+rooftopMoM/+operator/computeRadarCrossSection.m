function crossSection = computeRadarCrossSection(wavenumber, excitation, current, intensity)
% computeRadarCrossSection: Returns the monostatic radar cross section for the polarization of an incident plane
%   wave, (wavenumber * Z0)^2 * abs(excitation.' * current)^2 / (4 * pi * intensity^4).
%
%   Inputs:
%       wavenumber: Wavenumber of the wave.
%       excitation: Nx1 excitation vector of the plane wave.
%       current: Nx1 rooftop coefficients of the current induced by the wave.
%       intensity: Amplitude of the electric field of the wave in V/m.
%
%   Output:
%       crossSection: Radar cross section in m^2.
%
%   NOTE: The excitation integrates the rooftops against the incident wave, which is also the integral of the
%       backscattered far field, so the far field is proportional to excitation.' * current. Both of them scale with
%       the intensity, so the cross section does not depend on it.

    crossSection = (wavenumber * rooftopMoM.assembler.Constants.Z0)^2 * abs(excitation.' * current)^2 / ...
        (4 * pi * intensity^4);
end
