function [cX, cY] = generateGreensMatrixShape(primitive, direction)
% generateGreensMatrixShape: Calculate the number of rooftop elements in each direction for a given alignment.
%   This function determines the number of rooftop basis functions spanning the x and y directions.
%
%   Inputs:
%       primitive: Instance of the rooftopMoM.mesh.RectPrimitive class.
%
%       direction: Alignment of the rooftop basis functions.
%           'x' : Rooftops are oriented along the x-axis.
%           'y' : Rooftops are oriented along the y-axis.
%
%   Outputs:
%       cX: Number of rooftop elements spanning the x-directed edges.
%       cY: Number of rooftop elements spanning the y-directed edges.

    switch direction
        case 'x'
            % For x-directed rooftops:
            cX = primitive.countX() - 1;
            cY = primitive.countY();

        case 'y'
            % For y-directed rooftops:
            cX = primitive.countX();
            cY = primitive.countY() - 1;

        otherwise
            error("Unknown direction '%s'. Use 'x' or 'y'.", direction);
    end
end
