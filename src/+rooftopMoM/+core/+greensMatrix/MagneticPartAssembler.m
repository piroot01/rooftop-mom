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

        function chunkPlanner = chunkPlanner(obj)
            chunkPlanner = obj.assemblerChunkPlanner_;
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

            % TODO: Include the multi-wavenumber output
            G = obj.generateImplCore_(computeContext, chunks, slices, wavenumber)

            % Reshape according to the primitive shape.
            [cX, cY] = rooftopMoM.core.greensMatrix.utils.generateGreensMatrixShape(rooftops.primitive(), direction);

            % NOTE: the cY is aligned with columns and the cX with the rows
            G = reshape(G, cX, cY);
        end

        function G = generateImplCore_(~, context, chunks, slices, wavenumber)
            basisSize = context.basisSize;
            jacobians = context.jacobians;
            p1 = [1, 1, 2, 2];
            p2 = [1, 2, 1, 2];

            % Initialize greens matrix for this frequency.
            G = complex(zeros(basisSize, 1));

            % Process each chunk sequentially.
            for n = 1:numel(chunks)
                fprintf("Chunk ID=%d\n", n);

                chunk = chunks{n};
                slice = slices{n};

                % Determine row and column index ranges for current chunk.
                rowSlice = slice.rowIndicesStart:slice.rowIndicesEnd;
                colSlice = slice.columnIndicesStart:slice.columnIndicesEnd;

                rowRooftopPoints = context.rooftopPoints(:, rowSlice, :);
                columnRooftopPoints = context.rooftopPoints(:, colSlice, :);

                rowQuadraturePointsPolarity = context.quadraturePoints(:, rowSlice, p1);
                columnQuadraturePointsPolarity = context.quadraturePoints(:, colSlice, p2);

                distanceMatrix = rooftopMoM.core.impedanceMatrix.utils.generateDistanceMatrix(rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);
                integralKernel = rooftopMoM.core.impedanceMatrix.utils.generateScalarGreensFunction(wavenumber, distanceMatrix);

                % Obtain mask of singular interactions.
                singularityMask = context.singularityContext.indexMask(chunk.rowIndices, chunk.columnIndices, :);

                % Bundle chunk-specific data for integration.
                chunkContext = struct( ...
                    'chunkSize', [numel(chunk.rowIndices), numel(chunk.columnIndices)], ...
                    'rowRooftopPoints', rowRooftopPoints, ...
                    'columnRooftopPoints', columnRooftopPoints, ...
                    'singularityMask', singularityMask, ...
                    'distanceMatrix', distanceMatrix, ...
                    'integralKernel', integralKernel ...
                );

                % Compute non-singular and singular contributions for this chunk.
                [matrixChunk, singularityMaskChunkNonZeroDot] = rooftopMoM.core.impedanceMatrix.MagneticPartAssembler.generateChunk(context, chunkContext, wavenumber);
                singularCorrection = rooftopMoM.core.greensMatrix.MagneticPartAssembler.generateSingularities(context.singularityContext, chunkContext, singularityMaskChunkNonZeroDot);

                % Accumulate into the global Z matrix.
                G(chunk.rowIndices, chunk.columnIndices) = matrixChunk + singularCorrection;
            end

            % Apply final scaling by free-space impedance and Jacobian factors.
            G = 1i * wavenumber * rooftopMoM.core.Constants.Z0 / (4 * pi) .* G .* (jacobians .* jacobians);
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
 
        function singularities = generateSingularities(singularityContext, chunkContext, singularityMask)
            % Define polarity groupings for same and opposite polarities.
            samePolarities = [1, 4];
            differentPolarities = [2, 3];
            pCount = numel(samePolarities) + numel(differentPolarities);

            % Allocate polarity-resolved singular contributions.
            singularitiesPolarity = zeros([chunkContext.chunkSize, pCount]);
            singularitiesPolarity(:, :, samePolarities)      = singularityMask(:, :, samePolarities) .* singularityContext.type1Values;
            singularitiesPolarity(:, :, differentPolarities) = singularityMask(:, :, differentPolarities) .* singularityContext.type2Values;

            % Sum contributions across all polarities.
            singularities = sum(singularitiesPolarity, 3);
        end
    end
end
