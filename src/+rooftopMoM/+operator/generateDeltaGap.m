function excitation = generateDeltaGap(rooftops, port, voltage)
% generateDeltaGap: Returns the excitation vector of a delta-gap source, the voltage times the edge length on the
%   rooftops of the port and zero elsewhere.
%
%   Inputs:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%       port: 1xP indices of the rooftops across the gap.
%       voltage: Voltage of the source in V.
%
%   Output:
%       excitation: Nx1 vector.

    excitation = zeros(rooftops.count(), 1);
    excitation(port) = voltage * rooftops.edgeLengths(port);
end
