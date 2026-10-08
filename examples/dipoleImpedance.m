% dipoleImpedance: Input impedance of a strip dipole 1 m long and 4 mm wide, fed by a delta gap in its middle, over
%   frequency. The reactance vanishes a little below the half-wave length.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

dipoleLength = 1;
segmentCount = 40;
mesh = rooftopMoM.Mesh(segmentCount, 1, dipoleLength / segmentCount, 0.004);
solver = rooftopMoM.Solver('method', 'gmres', 'tolerance', 1e-8);

% Dipole lengths in wavelengths; every frequency is a problem of its own.
electricalLengths = 0.40:0.01:0.60;
impedances = zeros(size(electricalLengths));

for n = 1:numel(electricalLengths)
    problem = rooftopMoM.Problem(mesh, 2 * pi * electricalLengths(n) / dipoleLength);

    % The default port is the middle line of edges; 'port' takes any other rooftops.
    excitation = rooftopMoM.Excitation.deltaGap(problem, 'voltage', 1);
    impedances(n) = solver.solve(rooftopMoM.Operator(problem), excitation).inputImpedance();
end

plot(electricalLengths, real(impedances), electricalLengths, imag(impedances));
xlabel('dipole length / wavelength');
ylabel('input impedance (\Omega)');
legend('resistance', 'reactance');
grid on;

% The real part of the Green's matrix gives the radiated power, which equals the power fed into the gap.
problem = rooftopMoM.Problem(mesh, pi / dipoleLength);
result = solver.solve(rooftopMoM.Operator(problem), rooftopMoM.Excitation.deltaGap(problem));
current = result.solution();
resistance = rooftopMoM.Operator(problem, 'modifier', @real);

impedance = result.inputImpedance();
radiatedPower = real(current' * resistance.apply(current)) / 2;
fprintf('Half-wave dipole: %.1f %+.1fi ohm, radiates %.4f mW of the %.4f mW fed at 1 V.\n', real(impedance), ...
    imag(impedance), 1e3 * radiatedPower, 1e3 * real(1 / impedance) / 2);
