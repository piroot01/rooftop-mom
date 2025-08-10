classdef ElectricPartAssembler < rooftopMoM.core.greensMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_;
    end

    methods
        function obj = ElectricPartAssembler(assemblerChunkPlanner)
            obj.assemblerChunkPlanner_ = assemblerChunkPlanner;
        end
    end

    methods
        function G = generateFull(obj, context)
            % Produce all chunks for full matrix assembly
            chunks = obj.assemblerChunkPlanner_.generateFull(context)

            G = obj.generateSpecified(context, chunks);
        end

        function G = generateSpecified(obj, context, chunks)
            G = obj.generateImpl_(context, chunks);
        end
    end

    methods (Access = private)
        function G = generateImpl_(obj, context, chunks)
            wavenumber = context.wavenumber();

            rooftops = context.rooftops();
            quadrature = context.quadrature();

            basisSize = rooftops.countDir(dir);
            order = quadrature.order();
            blockSize = order^2;

            [points, values] = rooftops.generateRooftopsDir(dir);
            [quadraturePoints, rooftopPoints] = rooftopMoM.core.impedanceMatrix.utils.convertRooftopPointsToArrays(points);

            slices = rooftopMoM.core.impedanceMatrix.utils.generateSlicesFromChunks(chunks, order);

            jacobians = rooftopMoM.core.greensMatrix.utils.generate2DJacobians(rooftops.primitive(), basisSize);

            singularityContext = rooftopMoM.core.greensMatrix.MagneticPartAssembler.generateSingularityContext(rooftops.indicesDir(dir), values);

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

            % TODO: Include the multi-wavenumber output
            G = obj.generateImplCore_(computeContext, chunks, slices, wavenumber)

            % Reshape according to the primitive shape.
            [cX, cY] = rooftopMoM.core.greensMatrix.utils.generateGreensMatrixShape(rooftops.primitive(), dir);

            % NOTE: the cY is aligned with columns and the cX with the rows
            G = reshape(G, cX, cY);
        end
    end
end
