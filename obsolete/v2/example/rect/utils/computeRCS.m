function rho = computeRCS(wavenumber, current, excitation, waveIntensity)
    rho = (rooftopMoM.core.Constants.Z0 * wavenumber)^2 * abs(excitation.' * current)^2 / (4 * pi * waveIntensity^2);
end
