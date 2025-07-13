classdef UniformRooftops < rooftopMoM.basisFunctions.RectDivergenceConforming
    properties (Access = private)
        quadrature_ rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature;
    end

    methods
        function obj = UniformRooftops(mesh, quadrature)
            obj@rooftopMoM.basisFunctions.RectDivergenceConforming(mesh);
            obj.quadrature_ = quadrature;
        end
    end

    methods
        function generate(obj)
            generate@rooftopMoM.basisFunctions.RectDivergenceConforming(obj);

            obj.quadrature_.generate();
        end

        function points = generateQuadraturePoints(obj)
            % Retrieve basis function indices for all meshes.
            indices = obj.indicesCell(false);

            % Number of quadrature points per one basis function.
            pointCountPerBlock = obj.quadrature_.order()^2;

            % Total number of points per one page from two.
            pointCountPerPage = obj.count() * pointCountPerBlock;

            [points2Du, points2Dv] = obj.quadrature_.points();

            % Preallocate output.
            points = zeros(3, pointCountPerPage, 2, 'double');

            globalOffset = 1;
            meshCount = obj.mesh().count();

            for m = 1:meshCount
                % Obtain basis function indices for single mesh.
                localIndices = indices{m};

                basisSize = size(localIndices, 2);
                pointCountPerMesh = basisSize * pointCountPerBlock;

                % Extract single mesh geometry.
                vertices = obj.mesh().vertices{m};
                faces = obj.mesh().faces{m};

                % Preallocate local output.
                localPoints = zeros(3, pointCountPerMesh, 2, 'double');

                for n = 1:basisSize
                    % Determine slice indices for current block.
                    blockSlice = ((n - 1) * pointCountPerBlock + 1):(pointCountPerBlock * n);

                    % The rooftop faces share the same principal edges.
                    [e1, e2] = rooftopMoM.mesh.utils.generateRectFacePrincipalEdges(vertices, faces, localIndices(1, n));

                    % For each face of the rooftop compute its shift.
                    shifts = rooftopMoM.mesh.utils.generateRectFaceOrigin(vertices, faces, localIndices(:, n));

                    % Compute quadrature points for both rooftop faces simultaneously.
                    [pointsPositive, pointsNegative] = rooftopMoM.basisFunctions.utils.generateQuadraturePointsUniformRooftop(e1, e2, shifts, points2Du, points2Dv);

                    localPoints(:, blockSlice, 1) = pointsPositive;
                    localPoints(:, blockSlice, 2) = pointsNegative;
                end

                % Determine slice indices for current mesh.
                meshSlice = globalOffset:(globalOffset + pointCountPerMesh - 1);
                points(:, meshSlice, :) = localPoints;

                globalOffset = globalOffset + pointCountPerMesh;
            end
        end

        function points = generateRooftopPoints(obj)
        end

        function points = generateRooftopDivPoints(obj)
        end

        function quadrature = quadrature(obj)
            quadrature = obj.quadrature_;
        end
    end

    methods (Access = private)
    end
end
