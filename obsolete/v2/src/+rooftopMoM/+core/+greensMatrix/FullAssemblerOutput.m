classdef FullAssemblerOutput < handle
% FullAssemblerOutput: Stores and provides access to the computed Green's operator matrices
% for a given primitive in the form of a struct with fields: XX, XY, YX, YY

    properties (Access = private)
        primitive_;
        greensMatrix_;
    end

    methods
        function obj = FullAssemblerOutput(primitive)
            obj.primitive_ = primitive;
            obj.greensMatrix_ = struct('XX', [], 'XY', [], 'YX', [], 'YY', []);
        end

        function setGreensMatrix(obj, matrixStruct)
            requiredFields = {'XX', 'XY', 'YX', 'YY'};

            for f = requiredFields
                if ~isfield(matrixStruct, f{1})
                    error("The Green's matrix struct must contain field '%s'.", f{1});
                end
            end

            obj.greensMatrix_ = matrixStruct;
        end

        function setGreensMatrixField(obj, fieldName, matrix)
            % Sets a single Green's function matrix field
            obj.greensMatrix_.(fieldName) = matrix;
        end

        function matrix = getGreensMatrixField(obj, fieldName)
            % Returns the specified Green's function matrix field
            matrix = obj.greensMatrix_.(fieldName);
        end
    end

    methods
        function matrix = XX(obj), matrix = obj.greensMatrix_.XX; end
        function matrix = XY(obj), matrix = obj.greensMatrix_.XY; end
        function matrix = YX(obj), matrix = obj.greensMatrix_.YX; end
        function matrix = YY(obj), matrix = obj.greensMatrix_.YY; end
    end
end

