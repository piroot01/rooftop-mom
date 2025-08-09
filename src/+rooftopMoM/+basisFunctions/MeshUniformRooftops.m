classdef MeshUniformRooftops < rooftopMoM.basisFunctions.MeshDivergenceConforming
    properties (Access = private)
        quadrature_ rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature;
    end

    methods
        function obj = MeshUniformRooftops(mesh, quadrature)
            obj@rooftopMoM.basisFunctions.MeshDivergenceConforming(mesh);
            obj.quadrature_ = quadrature;
        end
    end

    methods
        function generate(obj)
            generate@rooftopMoM.basisFunctions.MeshDivergenceConforming(obj);
            obj.quadrature_.generate();
        end

        function [rooftopPoints, rooftopValues] = generateRooftops(obj)
            indices = obj.indices();
            meshCount = obj.mesh().count();
            rooftopPoints = cell(1, obj.count());
            rooftopValues = cell(1, obj.count());
            offset = uint32(0);

            for m = 1:meshCount
                localIndices = indices{m};
                basisSize = size(localIndices, 2);

                vertices = obj.mesh().vertices{m};
                faces = obj.mesh().faces{m};

                for n = 1:basisSize
                    [points, values] = obj.generateRooftop_(vertices, faces, localIndices(:, n));

                    k = offset + n;
                    rooftopPoints{k} = points;
                    rooftopValues{k} = values;
                end

                offset = offset + basisSize;
            end
        end

        function quadrature = quadrature(obj)
            quadrature = obj.quadrature_;
        end
    end

    methods (Access = private)
        function [rooftopPoints, rooftopValues] = generateRooftop_(obj, vertices, faces, facePair)
            [pointsU, pointsV] = obj.quadrature().points();

            [e1, e2] = rooftopMoM.mesh.utils.generateRectFacePrincipalEdges(vertices, faces, facePair(1));
            shifts = rooftopMoM.mesh.utils.generateRectFaceOrigin(vertices, faces, facePair);
            support = rooftopMoM.basisFunctions.utils.generateRooftopSupport(vertices, faces, facePair);

            [qpPos, qpNeg] = rooftopMoM.basisFunctions.utils.generateQuadraturePointsUniformRooftop(e1, e2, shifts, pointsU, pointsV);
            [rpPos, rpNeg] = rooftopMoM.basisFunctions.utils.generateRooftopPointsUniformRooftop(support, qpPos, qpNeg);
            [len, ht] = rooftopMoM.basisFunctions.utils.generateRooftopParametersUniformRooftop(support);

            % We separate the data about the rooftops into two struct due to the big memory requirements of the rooftopPoints.
            rooftopPoints = struct( ...
                'quadraturePointsPositive', qpPos, ...
                'quadraturePointsNegative', qpNeg, ...
                'rooftopPointsPositive', rpPos, ...
                'rooftopPointsNegative', rpNeg ...
            );

            rooftopValues = struct( ...
                'rooftopDivPositive', 1 / len, ...
                'rooftopDivNegative', (-1) / len, ...
                'length', len, ...
                'height', ht ...
            );
        end
    end
end
