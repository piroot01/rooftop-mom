classdef CirculantMatrixEmbedding < handle
    methods (Static)
        function result = compute(matrix, gridSize)
            M = gridSize(1);
            N = gridSize(2);

            % pad the input matrix to (M+1) x (N+1)
            matrixSize = size(matrix);
            coreMatrix = zeros(gridSize + 1);
            coreMatrix(1:matrixSize(1), 1:matrixSize(2)) = matrix;

            % precompute indices
            rowsTop = 1:M;
            rowsBottom = (M+1):(-1):2;
            colsLeft = 1:N;
            colsRight = (N+1):(-1):2;

            result = zeros(2 * M, 2 * N);

            % top-left quadrant
            result(1:M, 1:N) = coreMatrix(rowsTop, colsLeft);

            % bottom-left quadrant: reversed rows
            result(M+1:end, 1:N) = coreMatrix(rowsBottom, colsLeft);

            % top-right quadrant: reversed columns
            result(1:M, N+1:end) = coreMatrix(rowsTop, colsRight);

            % bottom-right quadrant: reversed rows and columns
            result(M+1:end, N+1:end) = coreMatrix(rowsBottom, colsRight);
        end
    end
end
