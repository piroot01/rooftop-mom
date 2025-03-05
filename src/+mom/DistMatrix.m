classdef DistMatrix
    methods (Static)
        function matrix = getChunk(quadPoints_i, quadPoints_j)
            % we use the fact that |a - b|^2 = |a|^2 + |b|^2 - 2a.b
            Si = sum(quadPoints_i.^2, 1);
            Sj = sum(quadPoints_j.^2, 1);

            dotProd = quadPoints_i.' * quadPoints_j;

            matrix = sqrt(Si.' + Sj - 2 * dotProd);
        end
    end
end