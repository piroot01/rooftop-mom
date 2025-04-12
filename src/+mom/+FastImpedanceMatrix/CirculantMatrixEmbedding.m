classdef CirculantMatrixEmbedding < handle
    methods (Static)
        function result = compute(matrix)
            gridSize = size(matrix) - 1;
            M = gridSize(1);
            N = gridSize(2);

            result = zeros(2 * M, 2 * N);

            % top-left block remains as in original matrix
            result(1:M, 1:N) = matrix(1:M, 1:N);

            % bottom-left block: reverse rows along the first dimension
            result((M + 1):end, 1:N) = matrix((M + 1):(-1):2, 1:N);

            % top-right block: reverse columns along the second dimension
            result(1:M, (N + 1):end) = matrix(1:M, (N + 1):(-1):2);

            % bottom-right block: reverse both rows and columns
            result((M + 1):end, (N + 1):end) = matrix((M + 1):(-1):2, (N + 1):(-1):2);
        end
    end
end
