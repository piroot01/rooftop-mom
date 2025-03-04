classdef MeshRect < handle
    properties (Access = private)
        rectPrimitive_ mom.MeshRectPrimitive;

        % origin must be a 3D vector
        origin_ double;
        rotX_ double = 0;
        rotY_ double = 0;
        rotZ_ double = 0;
    end

    methods
        function obj = MeshRect(rectPrimitive, origin, rotX, rotY, rotZ)
            obj.rectPrimitive_ = rectPrimitive;
            obj.origin_ = origin;
            obj.rotX_ = rotX;
            obj.rotY_ = rotY;
            obj.rotZ_ = rotZ;
        end
    end

    methods
        function mesh = getAsMesh(obj)
            function points = rectSide(count, size)
                rawPoints = 0:count;
                points = double(rawPoints) .* size;
            end

            function [vertices, transformation] = getVertices(pointsX, pointsY)
                [x, y] = meshgrid(pointsX, pointsY);
                points = [x(:), y(:), zeros(numel(x), 1)];

                transformationH = makehgtform('xrotate', obj.rotX_, 'yrotate', obj.rotY_, 'zrotate', obj.rotZ_);
                transformation = transformationH(1:3, 1:3);
                vertices = transformation * points.' + obj.origin_(:);
            end

            function faces = getFaces(countX, countY)
                [ii, jj] = meshgrid(1:countX, 1:countY);
                v1 = ii + (jj - 1) .* (countX + 1);
                v2 = ii + 1 + (jj - 1) .* (countX + 1);
                v3 = ii + 1 + jj .* (countX + 1);
                v4 = ii + jj * (countX + 1);
                faces = uint32([v1(:), v2(:), v3(:), v4(:)].');
            end

            countX = obj.rectPrimitive_.getPatchCountX();
            countY = obj.rectPrimitive_.getPatchCountY();
            sizeX = obj.rectPrimitive_.getPatchSizeX();
            sizeY = obj.rectPrimitive_.getPatchSizeY();

            pointsX = rectSide(countX, sizeX);
            pointsY = rectSide(countY, sizeY);

            mesh = mom.Mesh2D();
            mesh.origin = obj.origin_(:);
            [mesh.vertices, mesh.transformation] = getVertices(pointsX, pointsY);
            mesh.faces = getFaces(countX, countY);
       end
    end
end