function [quadraturePoints, rooftopPoints, rooftopValues] = generateRooftopData(rooftops, useCachedData)
% generateRooftopData: Returns arrays of quadrature points, rooftop sample points and their values.
%   Supports rooftopMoM.basisFunctions.MeshUniformRooftops and rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops.
%
%   Inputs:
%       rooftops: Instance of supported rooftops (MeshUniformRooftops or UniformRectPrimitiveUniformRooftops).
%       useCachedData (optional): logical scalar; if true, the function attempts to
%                 load precomputed quadrature/rooftop points and values from rooftopMoM.Cache.
%                 If required cached entries are missing an error is raised.
%
%   Outputs:
%       quadraturePoints: 3xNx2 array of quadrature point coordinates used for integration.
%       rooftopPoints:    3xNx2 (or 3xNx2 array of rooftop sample point coordinates corresponding to quadraturePoints.
%       rooftopValues:    1xN cell array (or 1x2 cell/struct array) containing per-rooftop auxiliary values returned by the rooftop generator.

    if nargin < 2
        useCachedData = false;
    end

    if useCachedData
        % Try to obtain the quadraturePoints.
        if rooftopMoM.Cache.getOptionDataStatus(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_QUADRATUREPOINTS)
            quadraturePoints = rooftopMoM.Cache.get(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_QUADRATUREPOINTS);
        elseif rooftopMoM.Cache.getOptionDataStatus(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_QUADRATUREPOINTS, 'X')) || ...
               rooftopMoM.Cache.getOptionDataStatus(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_QUADRATUREPOINTS, 'Y'))
            quadraturePoints = cat(2, ...
                rooftopMoM.Cache.get(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_QUADRATUREPOINTS, 'X')), ...
                rooftopMoM.Cache.get(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_QUADRATUREPOINTS, 'Y')) ...
            );
        else
            error('The quadraturePoints are not cached!');
        end

        % Try to obtain the rooftopPoints
        if rooftopMoM.Cache.getOptionDataStatus(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS)
            rooftopPoints = rooftopMoM.Cache.get(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS);
        elseif rooftopMoM.Cache.getOptionDataStatus(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPPOINTS, 'X')) || ...
               rooftopMoM.Cache.getOptionDataStatus(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPPOINTS, 'Y'))
            rooftopPoints = cat(2, ...
                rooftopMoM.Cache.get(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPPOINTS, 'X')), ...
                rooftopMoM.Cache.get(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPPOINTS, 'Y')) ...
            );
        else
            error('The rooftopPoints are not cached!');
        end

        % Try to obtain the rooftopValues
        if rooftopMoM.Cache.getOptionDataStatus(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPVALUES)
            rooftopValues = rooftopMoM.Cache.get(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPVALUES);
        elseif rooftopMoM.Cache.getOptionDataStatus(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPVALUES, 'X')) || ...
               rooftopMoM.Cache.getOptionDataStatus(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPVALUES, 'Y'))
            rooftopValues = cat(2, ...
                rooftopMoM.Cache.get(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPVALUES, 'X')), ...
                rooftopMoM.Cache.get(strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPVALUES, 'Y')) ...
            );
        else
            error('The rooftopValues are not cached!');
        end
    else
        if isa(rooftops, 'rooftopMoM.basisFunctions.MeshUniformRooftops')
            [points, rooftopValues] = rooftops.generateRooftops();
        elseif isa(rooftops, 'rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops')
            primitive = rooftops.primitive();
            [gridSizeX, gridSizeY] = primitive.gridSize();

            if gridSizeX ~= 1 && gridSizeY ~= 1
                [pointsX, valuesX] = rooftops.generateRooftopsDir('x');
                [pointsY, valuesY] = rooftops.generateRooftopsDir('y');
                points = cat(2, pointsX, pointsY);
                rooftopValues = cat(2, valuesX, valuesY);
            elseif gridSizeX ~= 1 && gridSizeY == 1
                [points, rooftopValues] = rooftops.generateRooftopsDir('x');
            elseif gridSizeX == 1 && gridSizeY ~= 1
                [points, rooftopValues] = rooftops.generateRooftopsDir('y');
            else
                % Both are 1 (degenerate) -> error (illegal primitive).
                error('The dimensions of the primitive are illegal %dx%d.', gridSizeX, gridSizeY);
            end
        end

        [quadraturePoints, rooftopPoints] = rooftopMoM.core.impedanceMatrix.utils.convertRooftopPointsToArrays(points);
    end
end
