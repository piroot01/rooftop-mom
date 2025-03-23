classdef DistanceMatrix
    methods (Static)
        function matrix = compute(rowQuadraturePoints, columnQuadraturePoints)
            % we use the fact that |a - b|^2 = |a|^2 + |b|^2 - 2a.b
            rowSum = sum(rowQuadraturePoints.^2, 1);
            columnSum = sum(columnQuadraturePoints.^2, 1);

            dotProduct = pagemtimes(permute(rowQuadraturePoints, [2, 1, 3]), columnQuadraturePoints);

            matrix = sqrt(permute(rowSum, [2, 1, 3]) + columnSum - 2 * dotProduct);
        end

        function matrix = getChunk(quadPoints_i, quadPoints_j)
            % we use the fact that |a - b|^2 = |a|^2 + |b|^2 - 2a.b
            Si = sum(quadPoints_i.^2, 1);
            Sj = sum(quadPoints_j.^2, 1);

            dotProd = quadPoints_i.' * quadPoints_j;

            matrix = sqrt(Si.' + Sj - 2 * dotProd);
        end
    end
end
