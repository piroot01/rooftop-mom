classdef UniformRect < handle
    properties (Access = private) 
        uniformRectPrimitive_ rooftopMoM.mesh.UniformRectPrimitive;

        origin_ double = [0; 0; 0];
        rotX_ double = 0;
        rotY_ double = 0;
        rotZ_ double = 0;
    end

    methods
        function obj = UniformRect(uniformRectPrimitive, origin, rotX, rotY, rotZ)
            obj.uniformRectPrimitive_ = uniformRectPrimitive;

            if nargin > 1
                obj.origin_ = origin(:);
                obj.rotX_   = rotX;
                obj.rotY_   = rotY;
                obj.rotZ_   = rotZ;
            end
        end
    end

    methods
        function mesh = generateMesh(obj, useGenome)
            if nargin < 2
                useGenome = true;
            end

            countX = obj.uniformRectPrimitive_.countX();
            countY = obj.uniformRectPrimitive_.countY();
            sizeX = obj.uniformRectPrimitive_.sizeX();
            sizeY = obj.uniformRectPrimitive_.sizeY();

            pointsX = obj.generatePointsPerSide_(countX, sizeX);
            pointsY = obj.generatePointsPerSide_(countY, sizeY);

            vertices = obj.generateVertices_(pointsX, pointsY);
            faces = obj.generateFaces_();

            if useGenome
                indicesGenome = obj.uniformRectPrimitive_.indicesGenome();
                facesGenome = faces(:, indicesGenome);

                % Identify which vertices are actually used by the culled faces.
                usedVertexIndices = unique(facesGenome(:));
                newIndices = zeros(1, size(vertices, 2), 'uint32');
                newIndices(usedVertexIndices) = uint32(1:numel(usedVertexIndices));

                % Remap face indices to the new compact numbering and remove unused vertices.
                vertices = vertices(:, usedVertexIndices);
                faces = newIndices(facesGenome);
            end

            mesh = rooftopMoM.mesh.Mesh(vertices, faces);
        end
    end

    methods (Access = private)
        function points = generatePointsPerSide_(~, count, spacing)
            indices = 0:count;
            points = double(indices) .* spacing;
        end

        function vertices = generateVertices_(obj, pointsX, pointsY)
            % Build a 2D grid of sample points in the local X–Y plane.
            [gridX, gridY] = meshgrid(pointsY, pointsX);
            localPoints = [gridX(:), gridY(:), zeros(numel(gridX), 1)];

            % Build the 3×3 rotation from the stored Euler angles.
            hgTransformation = makehgtform('xrotate', obj.rotX_, 'yrotate', obj.rotY_, 'zrotate', obj.rotZ_);
            rotationMatrix = hgTransformation(1:3, 1:3);

            % Apply rotation and then translate by the origin.
            vertices = rotationMatrix * localPoints.' + obj.origin_(:);
        end

        function faces = generateFaces_(obj)
            countX = obj.uniformRectPrimitive_.countX();
            countY = obj.uniformRectPrimitive_.countY();

            [gridX, gridY] = meshgrid(1:countX, 1:countY);

            % Compute the four corner‑vertex indices for each cell:
            %   v1: lower‑left, v2: lower‑right,
            %   v3: upper‑right, v4: upper‑left.
            v1 = (countX + 1) * (gridY.' - 1) + gridX.';
            v2 = v1 + 1;
            v3 = v1 + countX + 2;
            v4 = v1 + countX + 1;

            % Pack into a 4×N uint32 array, where each column is one face.
            faces = uint32([v1(:), v2(:), v3(:), v4(:)].');
        end
    end
end
