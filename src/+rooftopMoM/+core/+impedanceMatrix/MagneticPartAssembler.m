classdef MagneticPartAssembler < rooftopMoM.core.impedanceMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_;
    end

    methods
        function obj = MagneticPartAssembler(assemblerChunkPlanner)
            obj.assemblerChunkPlanner_ = assemblerChunkPlanner;
        end
    end

    methods
        function output = generateFull(obj, context)
            % Produce all chunks for full matrix assembly
            chunks = obj.assemblerChunkPlanner_.generateFull(context, true);

            output = obj.generateSpecified(context, chunks);
        end

        function output = generateSpecified(obj, context, chunks)
            wavenumber = context.wavenumber();
            matrixCount = numel(wavenumber);

            rooftops = context.rooftops();
            quadrature = context.quadrature();

            basisSize = rooftops.count();
            order = quadrature.order();
            blockSize = order^2;

            [points, values] = rooftops.generateRooftops();
            [quadraturePoints, rooftopPoints] = rooftopMoM.core.impedanceMatrix.utils.convertRooftopPointsToArrays(points);

            slices = rooftopMoM.core.impedanceMatrix.utils.generateSlicesFromChunks(chunks, order);

            ids = rooftopMoM.basisFunctions.utils.generateIDMask(rooftops);
            idMask = rooftopMoM.core.impedanceMatrix.utils.generateIDMask(ids);
            jacobians = rooftopMoM.core.impedanceMatrix.utils.generate2DJacobians(rooftops.mesh(), ids);

            singularityContext = rooftopMoM.core.impedanceMatrix.MagneticPartAssembler.generateSingularityContext(rooftops.indices(), values, basisSize);

            weights = quadrature.weights4D();

            computeContext = struct( ...
                'basisSize', basisSize, ...
                'blockSize', blockSize, ...
                'quadraturePoints', quadraturePoints, ...
                'rooftopPoints', rooftopPoints, ...
                'jacobians', jacobians, ...
                'idMask', idMask, ...
                'singularityContext', singularityContext, ...
                'weights', weights ...
            );

            % Pre-allocate full impedance tensor.
            Z = complex(zeros(basisSize, basisSize, matrixCount));

            % Assemble each frequency-dependent matrix.
            for n = 1:matrixCount
                Z(:, :, n) = obj.generateImpl_(computeContext, chunks, slices, wavenumber(n));
            end

            output = rooftopMoM.core.impedanceMatrix.AssemblerOutput(Z);
        end

        function chunkPlanner = chunkPlanner(obj)
            chunkPlanner = obj.assemblerChunkPlanner_;
        end
    end

    methods (Access = private)
        function Z = generateImpl_(~, context, chunks, slices, wavenumber)
            basisSize = context.basisSize;
            jacobians = context.jacobians;
            p1 = [1, 1, 2, 2];
            p2 = [1, 2, 1, 2];

            % Initialize impedance matrix for this frequency.
            Z = complex(zeros(basisSize, basisSize));

            % Process each chunk sequentially.
            for n = 1:numel(chunks)
                DEBUG("MagneticPartAssembler.generateImpl_: Assembling chunk %d/%d.", n, numel(chunks));

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

                % Extract ID and size masks for singularity correction.
                idMask = context.idMask(chunk.rowIndices, chunk.columnIndices);
                sizeMask = context.singularityContext.sizeMask(chunk.rowIndices, chunk.columnIndices);

                % Obtain mask of singular interactions.
                singularityMask = context.singularityContext.indexMask(chunk.rowIndices, chunk.columnIndices, :);

                % Bundle chunk-specific data for integration.
                chunkContext = struct( ...
                    'chunkSize', [numel(chunk.rowIndices), numel(chunk.columnIndices)], ...
                    'rowRooftopPoints', rowRooftopPoints, ...
                    'columnRooftopPoints', columnRooftopPoints, ...
                    'singularityMask', singularityMask, ...
                    'distanceMatrix', distanceMatrix, ...
                    'integralKernel', integralKernel, ...
                    'idMask', idMask, ...
                    'sizeMask', sizeMask ...
                );

                % Compute non-singular and singular contributions for this chunk.
                [matrixChunk, singularityMaskChunkNonZeroDot] = rooftopMoM.core.impedanceMatrix.MagneticPartAssembler.generateChunk(context, chunkContext, wavenumber);
                singularCorrection = rooftopMoM.core.impedanceMatrix.MagneticPartAssembler.generateSingularities(context.singularityContext, chunkContext, singularityMaskChunkNonZeroDot);

                % Accumulate into the global Z matrix.
                Z(chunk.rowIndices, chunk.columnIndices) = matrixChunk + singularCorrection;
            end

            % Apply final scaling by free-space impedance and Jacobian factors.
            Z = 1i * wavenumber * rooftopMoM.core.Constants.Z0 / (4 * pi) .* Z .* (jacobians.' * jacobians);
        end
    end

    methods (Static)
        function [chunk, singularityMaskChunkNonZeroDot] = generateChunk(context, chunkContext, wavenumber)
            % Replicate quadrature weights over each block.
            weights = repmat(context.weights, chunkContext.chunkSize);

            % Retrieve distance and kernel data.
            distanceMatrix = chunkContext.distanceMatrix;
            integralKernel = chunkContext.integralKernel;

            % Define polarity index maps.
            p1 = [1, 1, 2, 2];
            p2 = [1, 2, 1, 2];
            pCount = numel(p1);

            % Map rooftop points to each polarity combination.
            rowRooftopPointsPolarity    = chunkContext.rowRooftopPoints(:, :, p1);
            columnRooftopPointsPolarity = chunkContext.columnRooftopPoints(:, :, p2);

            % Create expanded mask identifying singular quadrature pairs.
            bigSingularityMask = repelem(chunkContext.singularityMask, context.blockSize, context.blockSize);

            % Identify zero-distance locations needing special treatment.
            distanceMatrixZero = (distanceMatrix < eps);

            % Compute dot products across polarity pairs.
            rooftopPointsDot = pagemtimes(permute(rowRooftopPointsPolarity, [2, 1, 3]), columnRooftopPointsPolarity);
            rooftopPointsDotZero = (rooftopPointsDot < eps);

            % Refine singularity mask to exclude perpendicular overlaps.
            bigSingularityMaskNonZeroDot = ~rooftopPointsDotZero & bigSingularityMask;
            singularityMaskChunkNonZeroDot = bigSingularityMaskNonZeroDot(1:context.blockSize:end, 1:context.blockSize:end, :);

            % Adjust kernel values at singular locations by subtracting 1/r term.
            integralKernel(bigSingularityMaskNonZeroDot) = integralKernel(bigSingularityMaskNonZeroDot) - 1 ./ distanceMatrix(bigSingularityMaskNonZeroDot);

            integrands = rooftopPointsDot .* integralKernel;

            % Substitute analytical limit for zero-distances..
            integrands(distanceMatrixZero) = rooftopPointsDot(distanceMatrixZero) * (-1i) * wavenumber;

            % Weight each integrand before summation.
            integrals = integrands .* weights;

            % Reshape to 5D for block summation over quadrature.
            integralsReshaped = reshape(integrals, context.blockSize, chunkContext.chunkSize(1), context.blockSize, chunkContext.chunkSize(2), pCount);
            integralsSummed = sum(sum(integralsReshaped, 1), 3);

            % Collapse polarity dimension to produce chunk matrix.
            chunkPolarity = reshape(integralsSummed, [chunkContext.chunkSize, pCount]);
            chunk = sum(chunkPolarity, 3);
        end

        function singularityContext = generateSingularityContext(indices, values, count)
            % Generate mask identifying singular interactions between rooftops.
            indexMask = rooftopMoM.core.impedanceMatrix.utils.generateSingularityMask(indices);

            % Compute unique rooftop size permutations and replicate over basis.
            [sizes, indicesPerm] = rooftopMoM.core.impedanceMatrix.utils.generateRooftopSizesUniquePerms(values);
            sizeMask = repmat(indicesPerm, 1, count);

            % Evaluate singular integral constants for each size case.
            type1Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type1(sizes(:, 1), sizes(:, 2));
            type2Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type2(sizes(:, 1), sizes(:, 2));

            singularityContext = struct( ...
                'indexMask', indexMask, ...
                'sizeMask', sizeMask, ...
                'type1Values', type1Values, ...
                'type2Values', type2Values ...
            );
        end

        function singularities = generateSingularities(singularityContext, chunkContext, singularityMask)
            % Define polarity groupings for same and opposite polarities.
            samePolarities = [1, 4];
            differentPolarities = [2, 3];
            pCount = numel(samePolarities) + numel(differentPolarities);

            % Distribute singular integral values according to rooftop sizes.
            type1Values = singularityContext.type1Values(chunkContext.sizeMask);
            type2Values = singularityContext.type2Values(chunkContext.sizeMask);

            % Identify valid self-interactions from ID mask.
            validIndices = (chunkContext.idMask > 0);
            combinedMask = bsxfun(@and, singularityMask, validIndices);

            % Allocate polarity-resolved singular contributions.
            singularitiesPolarity = zeros([chunkContext.chunkSize, pCount]);
            singularitiesPolarity(:, :, samePolarities)      = combinedMask(:, :, samePolarities) .* repmat(type1Values, [1, 1, 2]);
            singularitiesPolarity(:, :, differentPolarities) = combinedMask(:, :, differentPolarities) .* repmat(type2Values, [1, 1, 2]);

            % Sum contributions across all polarities.
            singularities = sum(singularitiesPolarity, 3);
        end
    end
end
