function memory = memoryForRooftopPoints(basisFunctions, quadrature)
% memoryForRooftopPoints: Returns the memory occupied by the rooftopPoints structure
%
%   Input:
%       basisFunctions: Instance of the rooftopMoM.basisFunctions.BasisFunctionsInterface
%       quadrature: Instance of the rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature
%
%   Output:
%       memory: Calculates the memory required for the rooftop points storage.

    % One rooftop double count, each point is 3d.
    singleRooftopCount = quadrature.order()^2 * 3;

    count = basisFunctions.count();

    % 2 for positive and negative rooftopPoints and 2 for positive and negative quadraturePoints
    memory = (4 * count * singleRooftopCount) * rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.Double);
end
