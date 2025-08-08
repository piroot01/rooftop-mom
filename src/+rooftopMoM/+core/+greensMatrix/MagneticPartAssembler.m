classdef MagneticPartAssembler < rooftopMoM.core.greensMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_;
    end

    methods
        function obj = MagneticPartAssembler(assemblerChunkPlanner)
            obj.assemblerChunkPlanner_ = assemblerChunkPlanner;
        end
    end

    methods
        function [Gxx, Gyy] = generateFull(obj, context)
            % Produce all chunks for full matrix assembly
            chunksX = obj.assemblerChunkPlanner_.generateFullDirection(context, 'x');
            chunksY = obj.assemblerChunkPlanner_.generateFullDirection(context, 'y');

            [Gxx, Gyy] = obj.generateSpecified(context, chunksX, chunksY);
        end

        function [Gxx, Gyy] = generateSpecified(obj, context, chunksX, chunksY)
            Gxx = obj.generateDirectionSpecified(context, chunksX, 'x');
            Gyy = obj.generateDirectionSpecified(context, chunksY, 'y');
        end

        function G = generateFullDirection(obj, context, direction)
            % Produce all chunks for full matrix assembly
            chunks = obj.assemblerChunkPlanner_.generateFullDirection(context, direction);

            G = obj.generateImpl_(context, chunks, direction);
        end

        function G = generateSpecifiedDirection(obj, context, chunks, direction)
            G = obj.generateImpl_(context, chunks, direction);
        end
    end

    methods (Access = private)
        function G = generateImpl_(obj, context, chunks, direction)
            wavenumber = context.wavenumber();

            rooftops = context.rooftops();
            quadrature = context.quadrature();

            basisSize = rooftops.countDirection(direction);
            order = quadrature.order();
            blockSize = order^2;

            [points, values] = rooftops.generateRooftopsDirection(direction);
            [quadraturePoints, rooftopPoints] = rooftopMoM.core.impedanceMatrix.utils.convertRooftopPointsToArrays(points);

            slices = rooftopMoM.core.impedanceMatrix.utils.generateSlicesFromChunks(chunks, order);

            jacobians = rooftopMoM.core.greensMatrix.utils.generate2DJacobians(rooftops.primitive(), basisSize);

            singularityContext = rooftopMoM.core.greensMatrix.MagneticPartAssembler.generateSingularityContext(rooftops.indicesDirection(direction), values);

            weights = quadrature.weights4D();

            computeContext = struct( ...
                'basisSize', basisSize, ...
                'blockSize', blockSize, ...
                'quadraturePoints', quadraturePoints, ...
                'rooftopPoints', rooftopPoints, ...
                'jacobians', jacobians, ...
                'singularityContext', singularityContext, ...
                'weights', weights ...
            );

            % Pre-allocate full greens matrix.
            G = complex(zeros(basisSize, 1));

            % TODO: Call the generateImplCore_
        end

        function G = generateImplCore_(obj, context, chunks, slices, wavenumber)
        end
    end

    methods (Static)
        function singularityContext = generateSingularityContext(indices, values)
            % Generate mask identifying singular interactions between rooftops.
            indexMask = rooftopMoM.core.greensMatrix.utils.generateSingularityMask(indices);

            length = values.length;
            height = values.height;

            % Evaluate singular integral constants for each size case.
            type1Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type1(length, height);
            type2Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type2(length, height);

            singularityContext = struct( ...
                'indexMask', indexMask, ...
                'type1Values', type1Values, ...
                'type2Values', type2Values ...
            );
        end
    end
end
