classdef UniformRooftops
% UniformRooftops: Rooftop basis functions of unit height on a uniform rectangular grid.
%   Each interior edge carries one rooftop. It rises linearly from 0 to 1 across the face before the edge and falls
%   back to 0 across the face after it. The x-directed rooftops come first, then the y-directed ones; within one
%   direction the x index runs fastest.
%
%   NOTE: All rooftops of one direction are translates of a reference rooftop at the origin. This is what makes the
%       Green's matrix block Toeplitz.

    properties (Access = private)
        primitive_ rooftopMoM.mesh.UniformRectPrimitive;
        quadrature_ rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature;
    end

    methods
        function obj = UniformRooftops(primitive, quadrature)
            obj.primitive_ = primitive;
            obj.quadrature_ = quadrature;
        end
    end

    methods
        function primitive = primitive(obj)
            primitive = obj.primitive_;
        end

        function quadrature = quadrature(obj)
            quadrature = obj.quadrature_;
        end

        function count = count(obj)
        % count: Returns the number of all rooftops.

            count = obj.countDir('x') + obj.countDir('y');
        end

        function count = countDir(obj, direction)
        % countDir: Returns the number of rooftops of one direction, 'x' or 'y'.

            [rowCount, columnCount] = obj.gridSize(direction);
            count = rowCount * columnCount;
        end

        function [rowCount, columnCount] = gridSize(obj, direction)
        % gridSize: Returns the grid of the rooftops of one direction, (countX - 1) x countY for 'x' and
        %   countX x (countY - 1) for 'y'.

            if direction == 'x'
                rowCount = obj.primitive_.countX() - 1;
                columnCount = obj.primitive_.countY();
            else
                rowCount = obj.primitive_.countX();
                columnCount = obj.primitive_.countY() - 1;
            end
        end

        function content = content(obj)
        % content: Returns the RooftopContent of the grid.

            hasX = obj.countDir('x') > 0;
            hasY = obj.countDir('y') > 0;

            if hasX && hasY
                content = rooftopMoM.basisFunctions.RooftopContent.BOTH;
            elseif hasX
                content = rooftopMoM.basisFunctions.RooftopContent.XONLY;
            elseif hasY
                content = rooftopMoM.basisFunctions.RooftopContent.YONLY;
            else
                content = rooftopMoM.basisFunctions.RooftopContent.NONE;
            end
        end

        function offsets = faceOffsets(obj, direction, indices)
        % faceOffsets: Returns the position of the first face of the given rooftops in face sizes.
        %
        %   Inputs:
        %       direction: 'x' or 'y'.
        %       indices: 1xN vector of the rooftop indices within the direction.
        %
        %   Output:
        %       offsets: 2xN matrix of the face offsets (i; j), counted from 0.

            [rowCount, columnCount] = obj.gridSize(direction);
            [i, j] = ind2sub([rowCount, columnCount], indices);
            offsets = [i; j] - 1;
        end

        function rooftop = referenceRooftop(obj, direction)
        % referenceRooftop: Returns the rooftop of one direction whose first face is at the origin.
        %
        %   Input:
        %       direction: 'x' or 'y'.
        %
        %   Output:
        %       rooftop: Struct with the fields:
        %           secondFaceOffset: 2x1 offset of the second face from the first one.
        %           length: Face size along the direction.
        %           height: Face size across the direction, the length of the edge.
        %           amplitudes: Bx2 values at the B quadrature points of the first (rising) and the second
        %               (falling) face.

            [pointsU, pointsV] = obj.quadrature_.points();

            if direction == 'x'
                coordinate = pointsU.';
                rooftop.secondFaceOffset = [1; 0];
                rooftop.length = obj.primitive_.sizeX();
                rooftop.height = obj.primitive_.sizeY();
            else
                coordinate = pointsV.';
                rooftop.secondFaceOffset = [0; 1];
                rooftop.length = obj.primitive_.sizeY();
                rooftop.height = obj.primitive_.sizeX();
            end

            rooftop.amplitudes = [coordinate, 1 - coordinate];
        end
    end
end
