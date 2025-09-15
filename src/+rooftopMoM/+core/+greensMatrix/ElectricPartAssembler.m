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
            chunks = obj.assemblerChunkPlanner_.generateFull(context);

            G = obj.generateSpecified(context, chunks);
        end

        function G = generateSpecified(obj, context, chunks)
            % NOTE: Just to be consistent.
            G = obj.generateImpl_(context, chunks);
        end
    end

    methods (Access = private)
        function G = generateImpl_(obj, context, chunks)
            wavenumber = context.wavenumber();

            rooftops = context.rooftops();
            primitive = rooftops.primitive();
            quadrature = context.quadrature();

            count = primitive.count();
            order = quadrature.order();
            blockSize = order^2;

            quadraturePoints = rooftopMoM.basisFunctions.utils.generateQuadraturePointsUniformRectPrimitive(quadrature, primitive);

            slices = rooftopMoM.core.impedanceMatrix.utils.generateSlicesFromChunks(chunks, order);

            jacobians = rooftopMoM.core.greensMatrix.utils.generate2DJacobians(primitive, count);

            singularityContext = rooftopMoM.core.greensMatrix.ElectricPartAssembler.generateSingularityContext(primitive);

            weights = quadrature.weights4D();

            computeContext = struct( ...
                'count', count, ...
                'blockSize', blockSize, ...
                'quadraturePoints', quadraturePoints, ...
                'jacobians', jacobians, ...
                'singularityContext', singularityContext, ...
                'weights', weights ...
            );

            % TODO: Include the multi-wavenumber output
            G = obj.generateImplCore_(computeContext, chunks, slices, wavenumber);

            % Reshape according to the primitive shape.
            [cX, cY] = primitive.gridSize();

            % NOTE: the cY is aligned with columns and the cX with the rows
            G = reshape(G, cX, cY);
        end

        function G = generateImplCore_(~, context, chunks, slices, wavenumber)
            count = context.count;
            jacobians = context.jacobians;

            % Initialize greens matrix for this frequency.
            G = complex(zeros(count, 1));

            % Process each chunk sequentially.
            for n = 1:numel(chunks)
                DEBUG("ElectricPartAssembler.generateImplCore_: Assembling chunk %d/%d.", n, numel(chunks));

                chunk = chunks{n};
                slice = slices{n};

                % Determine row and column index ranges for current chunk.
                rowSlice = slice.rowIndicesStart:slice.rowIndicesEnd;
                colSlice = slice.columnIndicesStart:slice.columnIndicesEnd;

                rowQuadraturePointsPolarity = context.quadraturePoints(:, rowSlice);
                columnQuadraturePointsPolarity = context.quadraturePoints(:, colSlice);

                distanceMatrix = rooftopMoM.core.impedanceMatrix.utils.generateDistanceMatrix(rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);
                integralKernel = rooftopMoM.core.impedanceMatrix.utils.generateScalarGreensFunction(wavenumber, distanceMatrix);

                singularityMask = context.singularityContext.indexMask(chunk.rowIndices, chunk.columnIndices, :);

                % Bundle chunk-specific data for integration.
                chunkContext = struct( ...
                    'chunkSize', [numel(chunk.rowIndices), numel(chunk.columnIndices)], ...
                    'singularityMask', singularityMask, ...
                    'distanceMatrix', distanceMatrix, ...
                    'integralKernel', integralKernel ...
                );

                G(chunk.rowIndices, chunk.columnIndices) = rooftopMoM.core.greensMatrix.ElectricPartAssembler.generateChunk(context, chunkContext, wavenumber);
            end

            % Apply the singularity correction.
            G(context.singularityContext.indexMask) = G(context.singularityContext.indexMask) + context.singularityContext.type3Values;

            % Apply final scaling by free-space impedance and Jacobian factors.
            G = 1i * rooftopMoM.core.Constants.Z0 / (4 * pi * wavenumber) .* G .* (jacobians .* jacobians);
        end
    end

    methods (Static)
        function chunk = generateChunk(context, chunkContext, wavenumber)
            % Replicate quadrature weights over each block.
            weights = repmat(context.weights, chunkContext.chunkSize);

            % Retrieve distance and kernel data.
            distanceMatrix = chunkContext.distanceMatrix;
            integralKernel = chunkContext.integralKernel;

            % Create expanded mask identifying singular quadrature pairs.
            bigSingularityMask = repelem(chunkContext.singularityMask, context.blockSize, context.blockSize);

            % Identify zero-distance locations needing special treatment.
            distanceMatrixZero = (distanceMatrix < eps);

            % Adjust kernel values at singular locations by subtracting 1/r term.
            integralKernel(bigSingularityMask) = integralKernel(bigSingularityMask) - 1 ./ distanceMatrix(bigSingularityMask);

            % Substitute analytical limit for zero-distances..
            integralKernel(distanceMatrixZero) = (-1i) * wavenumber;

            % Weight each integrand before summation.
            integrals = integralKernel .* weights;

            % Reshape to 5D for block summation over quadrature.
            integralsReshaped = reshape(integrals, context.blockSize, chunkContext.chunkSize(1), context.blockSize, chunkContext.chunkSize(2));
            integralsSummed = sum(sum(integralsReshaped, 1), 3);

            % Collapse polarity dimension to produce chunk matrix.
            chunk = reshape(integralsSummed, [chunkContext.chunkSize]);
        end

        function singularityContext = generateSingularityContext(primitive)
            % NOTE: The singular value is 1 and located at the (1, 1) position of the electric part of the Green's matrix.
            indexMask = false(primitive.count(), 1);
            indexMask(1) = true;

            length = primitive.sizeX();
            height = primitive.sizeY();

            % Evaluate the singular integral constants for the given size.
            type3Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type3(length, height);

            singularityContext = struct( ...
                'indexMask', indexMask, ...
                'type3Values', type3Values ...
            );
        end
    end
end
