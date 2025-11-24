function lengths = generateRooftopEdgeLengths(rooftops, rooftopValues)
% generateRooftopEdgeLengths: Computes per-rooftop edge lengths for supported rooftop classes.
%
%   Inputs:
%       rooftops: Instance of MeshUniformRooftops or UniformRectPrimitiveUniformRooftops.
%       rooftopValues: Cell/struct array returned by generateRooftopData(), containing
%                      per-direction rooftop metadata including the edge length.
%
%   Outputs:
%       lengths: 1xN array of rooftop edge lengths. For MeshUniformRooftops lengths are extracted
%                from each rooftopValue entry. For UniformRectPrimitiveUniformRooftops the X and Y
%                lengths are expanded according to the number of rooftops in each direction.

    if isa(rooftops, 'rooftopMoM.basisFunctions.MeshUniformRooftops')
        lengths = cellfun(@(c) c.length, rooftopValues, 'UniformOutput', true);
    elseif isa(rooftops, 'rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops')
        basisSizeX = rooftops.countDir('x');
        basisSizeY = rooftops.countDir('y');
        lengths = cat(2, repelem(rooftopValues(1).length, 1, basisSizeX), repelem(rooftopValues(2).length, 1, basisSizeY));
    else
        error('Unsupported rooftop kind.');
    end
end
