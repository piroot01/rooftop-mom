classdef UniformRectPrimitiveDivergenceConforming < rooftopMoM.basisFunctions.UniformRectPrimitiveBasisFunctionsInterface
    properties (Access = private)
        primitive_;

        indicesX_;
        indicesY_;
    end

    methods
        function obj = UniformRectPrimitiveDivergenceConforming(primitive)
            obj.primitive_ = primitive;
        end
    end

    methods
        function generate(obj)
            [countX, countY] = obj.primitive_.gridSize();

            obj.indicesX_ = rooftopMoM.mesh.utils.generateAllInteriorEdgesInDir(countX, countY, 'x');
            obj.indicesY_ = rooftopMoM.mesh.utils.generateAllInteriorEdgesInDir(countX, countY, 'y');
        end

        function count = count(obj)
            count = obj.countDir('x') + obj.countDir('y');
        end

        function count = countDir(obj, dir)
            switch dir
                case 'x'
                    count = size(obj.indicesX_, 2);
                case 'y'
                    count = size(obj.indicesY_, 2);
                otherwise
                    error("Unknown direction '%s'. Use 'x' or 'y'.", direction);
            end
        end

        function indices = indices(obj)
            indices = cat(2, obj.indicesX_, obj.indicesY_);
        end

        function indices = indicesDir(obj, dir)
            switch dir
                case 'x'
                    indices = obj.indicesX_;
                case 'y'
                    indices = obj.indicesY_;
                otherwise
                    error("Unknown direction '%s'. Use 'x' or 'y'.", direction);
            end
        end

        function primitive = primitive(obj)
            primitive = obj.primitive_;
        end

        function memory = memory(obj)
            % TODO: Implement the memory requirements of this class.
            memory = nan;
        end
    end
end
