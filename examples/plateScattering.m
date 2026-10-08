% plateScattering: Monostatic radar cross section of a square plate two wavelengths wide, as the incident plane wave
%   tilts from normal incidence towards the x axis. Physical optics, a good estimate near normal incidence, is
%   plotted for comparison.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
rooftopMoM.Logging.configure('level', 'info');

wavelength = 1;
side = 2 * wavelength;
faceCount = 20;
problem = rooftopMoM.Problem(rooftopMoM.Mesh(faceCount, faceCount, side / faceCount, side / faceCount), ...
    2 * pi / wavelength);

% The product and the preconditioner are built once and serve every angle.
operator = rooftopMoM.Operator(problem);
preconditioner = rooftopMoM.Preconditioner(problem);
solver = rooftopMoM.Solver('tolerance', 1e-6);

angles = 0:2:40;
crossSections = zeros(size(angles));
current = [];

for n = 1:numel(angles)
    % The electric field stays along y, parallel to the plate.
    excitation = rooftopMoM.Excitation.planeWave(problem, 'direction', [sind(angles(n)); 0; -cosd(angles(n))], ...
        'polarization', [0; 1; 0]);

    % The current of the previous angle is a close initial guess.
    result = solver.solve(operator, excitation, 'preconditioner', preconditioner, 'initialGuess', current);
    current = result.solution();
    crossSections(n) = result.radarCrossSection();
end

rooftopMoM.Logging.reset();

% Physical optics: 4 * pi * (area / wavelength)^2 * cos^2, times sinc^2 of the phase across the plate.
phase = 2 * pi / wavelength * side * sind(angles);
sincs = sin(phase) ./ phase;
sincs(phase == 0) = 1;
opticsCrossSections = 4 * pi * (side^2 / wavelength)^2 * cosd(angles).^2 .* sincs.^2;

plot(angles, 10 * log10(crossSections), 'o', angles, 10 * log10(opticsCrossSections), '-');
ylim([-20, 30]);
xlabel('angle of incidence (deg)');
ylabel('radar cross section (dBsm)');
legend('rooftopMoM', 'physical optics');
grid on;
