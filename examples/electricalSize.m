% electricalSize: Radar cross section of square plates at normal incidence over their electrical size ka, with a the
%   radius of the smallest sphere around the plate. Every plate is meshed with 10 faces per wavelength, so each face
%   has the same electrical size at every ka, and sigma / lambda^2 compares the plates whatever their physical size.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

electricalSizes = 1:0.5:10;
crossSections = zeros(size(electricalSizes));
solver = rooftopMoM.Solver('tolerance', 1e-6);

for n = 1:numel(electricalSizes)
    problem = rooftopMoM.Problem.fromElectricalSize(electricalSizes(n), 'density', 10);
    result = solver.solve(rooftopMoM.Operator(problem), rooftopMoM.Excitation.planeWave(problem), ...
        'preconditioner', rooftopMoM.Preconditioner(problem));
    crossSections(n) = result.normalizedRadarCrossSection();
end

% Physical optics at normal incidence, 4 * pi * (area / lambda^2)^2, with the side of the square ka / (pi * sqrt(2))
% wavelengths.
opticsCrossSections = 4 * pi * (electricalSizes.^2 / (2 * pi^2)).^2;

plot(electricalSizes, 10 * log10(crossSections), 'o', electricalSizes, 10 * log10(opticsCrossSections), '-');
xlabel('electrical size ka');
ylabel('radar cross section / \lambda^2 (dB)');
legend('rooftopMoM', 'physical optics', 'Location', 'southeast');
grid on;
