classdef Inversion < handle
    methods (Static)
        function result = solve(impedanceMatrix, excitationVector)
            % For a single solve, this is the most recommended and efficient method.
            % MATLAB's `\` (mldivide) automatically detects that the matrix is
            % square and uses a highly optimized LU decomposition internally.
            result = impedanceMatrix \ excitationVector(:);
        end
    end
end
