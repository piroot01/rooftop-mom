classdef PerturbedMatrix
    properties (Access = private)
        unperturbedMatrix_;
        perturbedMatrix_;

        indices_ uint32;
    end

    methods
        function obj = PerturbedMatrix(matrix)
            obj.unperturbedMatrix_ = matrix;
            colCount = size(matrix, 2);
            obj.indices_ = 1:colCount;
        end
    end

    methods
        function matrix = unperturbedMatrix(obj)
            matrix = obj.unperturbedMatrix_;
        end

        function matrix = perturbedMatrix(obj)
            matrix = obj.perturbedMatrix_;
        end

        function indices = indices(obj)
            indices = obj.indices_;
        end
    end
end
