classdef Mesh
% Mesh: Uniform grid of countX x countY rectangular faces of the size sizeX x sizeY in m. It lies in the plane z = 0
%   with one corner at the origin, and every interior edge carries a rooftop basis function.
%       mesh = rooftopMoM.Mesh(64, 32, 0.05, 0.05);

    properties (Access = private)
        primitive_ rooftopMoM.mesh.UniformRectPrimitive;
    end

    methods
        function obj = Mesh(countX, countY, sizeX, sizeY)
        % Mesh: Creates the grid.
        %
        %   Inputs:
        %       countX: Number of faces in the x direction.
        %       countY: Number of faces in the y direction.
        %       sizeX: Face size in the x direction in m.
        %       sizeY: Face size in the y direction in m.

            arguments
                countX (1, 1) {mustBeInteger, mustBePositive}
                countY (1, 1) {mustBeInteger, mustBePositive}
                sizeX (1, 1) {mustBeReal, mustBeFinite, mustBePositive}
                sizeY (1, 1) {mustBeReal, mustBeFinite, mustBePositive}
            end

            if countX == 1 && countY == 1
                error('rooftopMoM:mesh:noRooftops', 'A single face has no interior edge, so it carries no rooftop.');
            end

            obj.primitive_ = rooftopMoM.mesh.UniformRectPrimitive(double(countX), double(countY), double(sizeX), ...
                double(sizeY));
        end
    end

    methods
        function [countX, countY] = count(obj)
        % count: Returns the numbers of faces in the x and the y direction.

            countX = obj.primitive_.countX();
            countY = obj.primitive_.countY();
        end

        function [sizeX, sizeY] = faceSize(obj)
        % faceSize: Returns the face sizes in the x and the y direction in m.

            sizeX = obj.primitive_.sizeX();
            sizeY = obj.primitive_.sizeY();
        end

        function radius = radius(obj)
        % radius: Returns the radius a of the smallest sphere around the plate, half its diagonal, in m. A problem on
        %   the mesh has the electrical size ka.

            radius = hypot(obj.primitive_.countX() * obj.primitive_.sizeX(), ...
                obj.primitive_.countY() * obj.primitive_.sizeY()) / 2;
        end

        function primitive = primitive(obj)
            primitive = obj.primitive_;
        end
    end
end
