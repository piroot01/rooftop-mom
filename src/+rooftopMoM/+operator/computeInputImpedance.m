function impedance = computeInputImpedance(rooftops, current, port, voltage)
% computeInputImpedance: Returns the input impedance of a delta-gap source, its voltage over the total current
%   through the edges of the port.
%
%   Inputs:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%       current: Nx1 rooftop coefficients of the current induced by the source.
%       port: 1xP indices of the rooftops across the gap.
%       voltage: Voltage of the source in V.
%
%   Output:
%       impedance: Input impedance in ohms.
%
%   NOTE: The rooftops have unit height, so a coefficient is the current density in A/m across the edge. Times the
%       edge length it gives the current through the edge.

    impedance = voltage / (rooftops.edgeLengths(port) * current(port));
end
