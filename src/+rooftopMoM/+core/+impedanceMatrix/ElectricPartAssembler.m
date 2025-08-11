classdef ElectricPartAssembler < rooftopMoM.core.impedanceMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_;
    end

    methods
        function obj = ElectricPartAssembler(assemblerChunkPlanner)
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
            quadraturePoints = rooftopMoM.core.impedanceMatrix.utils.stackPositiveAndNegativeFields(points, 'quadraturePointsPositive', 'quadraturePointsNegative');
            rooftopDivPoints = rooftopMoM.core.impedanceMatrix.utils.stackPositiveAndNegativeFields(values, 'rooftopDivPositive', 'rooftopDivNegative');

            slices = rooftopMoM.core.impedanceMatrix.utils.generateSlicesFromChunks(chunks, order);

            ids = rooftopMoM.basisFunctions.utils.generateIDMask(rooftops);
            idMask = rooftopMoM.core.impedanceMatrix.utils.generateIDMask(ids);
            jacobians = rooftopMoM.core.impedanceMatrix.utils.generate2DJacobians(rooftops.mesh(), ids);

            singularityContext = rooftopMoM.core.impedanceMatrix.ElectricPartAssembler.generateSingularityContext(rooftops.indices(), values, basisSize);

            singularities = rooftopMoM.core.impedanceMatrix.ElectricPartAssembler.generateSingularities(singularityContext, idMask);

            weights = quadrature.weights4D();

            computeContext = struct( ...
                'basisSize', basisSize, ...
                'blockSize', blockSize, ...
                'quadraturePoints', quadraturePoints, ...
                'rooftopDivPoints', rooftopDivPoints, ...
                'jacobians', jacobians, ...
                'idMask', idMask, ...
                'singularityContext', singularityContext, ...
                'singularities', singularities, ...
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
                fprintf("Chunk ID=%d\n", n);

                chunk = chunks{n};
                slice = slices{n};

                % Determine row and column index ranges for current chunk.
                rowSlice = slice.rowIndicesStart:slice.rowIndicesEnd;
                colSlice = slice.columnIndicesStart:slice.columnIndicesEnd;

                rowRooftopDivPoints = context.rooftopDivPoints(:, chunk.rowIndices, :);
                columnRooftopDivPoints = context.rooftopDivPoints(:, chunk.columnIndices, :);

                rowQuadraturePointsPolarity = context.quadraturePoints(:, rowSlice, p1);
                columnQuadraturePointsPolarity = context.quadraturePoints(:, colSlice, p2);

                distanceMatrix = rooftopMoM.core.impedanceMatrix.utils.generateDistanceMatrix(rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);
                integralKernel = rooftopMoM.core.impedanceMatrix.utils.generateScalarGreensFunction(wavenumber, distanceMatrix);

                % Obtain mask and values of singular interactions.
                singularityMask = context.singularityContext.indexMask(chunk.rowIndices, chunk.columnIndices, :);
                singularities = context.singularities(chunk.rowIndices, chunk.columnIndices, :);

                % Bundle chunk-specific data for integration.
                chunkContext = struct( ...
                    'chunkSize', [numel(chunk.rowIndices), numel(chunk.columnIndices)], ...
                    'rowRooftopDivPoints', rowRooftopDivPoints, ...
                    'columnRooftopDivPoints', columnRooftopDivPoints, ...
                    'singularityMask', singularityMask, ...
                    'electricSingularities', singularities, ...
                    'distanceMatrix', distanceMatrix, ...
                    'integralKernel', integralKernel ...
                );

                % Compute contributions for this chunk.
                Z(chunk.rowIndices, chunk.columnIndices) = rooftopMoM.core.impedanceMatrix.ElectricPartAssembler.generateChunk(context, chunkContext, wavenumber);
            end

            % Apply final scaling by free-space impedance and Jacobian factors.
            Z = 1i * rooftopMoM.core.Constants.Z0 / (4 * pi * wavenumber) .* Z .* (jacobians.' * jacobians);
        end
    end

    methods (Static)
        function chunk = generateChunk(context, chunkContext, wavenumber)
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
            rowRooftopDivPointsPolarity = chunkContext.rowRooftopDivPoints(:, :, p1);
            columnRooftopDivPointsPolarity = chunkContext.columnRooftopDivPoints(:, :, p2);

            % Create expanded mask identifying singular quadrature pairs.
            bigSingularityMask = repelem(chunkContext.singularityMask, context.blockSize, context.blockSize);

            % Identify zero-distance locations needing special treatment.
            distanceMatrixZero = (distanceMatrix < eps);

            % Compute dot products across polarity pairs.
            rooftopDivDot = pagemtimes(permute(rowRooftopDivPointsPolarity, [2, 1, 3]), columnRooftopDivPointsPolarity);

            % Adjust kernel values at singular locations by subtracting 1/r term.
            integralKernel(bigSingularityMask) = integralKernel(bigSingularityMask) - 1 ./ distanceMatrix(bigSingularityMask);

            % Substitute analytical limit for zero-distances..
            integralKernel(distanceMatrixZero) = (-1i) * wavenumber;

            % Weight each integrand before summation.
            integrals = integralKernel .* weights;

            % Reshape to 5D for block summation over quadrature.
            integralsReshaped = reshape(integrals, context.blockSize, chunkContext.chunkSize(1), context.blockSize, chunkContext.chunkSize(2), pCount);
            integralsSummed = sum(sum(integralsReshaped, 1), 3);

            % Collapse polarity dimension to produce chunk matrix.
            chunkPolarity = reshape(integralsSummed, [chunkContext.chunkSize, pCount]);

            % Apply the singular contributions.
            chunkPolarity = chunkPolarity + chunkContext.electricSingularities;

            % Scale by the rooftop divergences.
            chunkPolarity = chunkPolarity .* rooftopDivDot;

            chunk = sum(chunkPolarity, 3);
        end

        function singularityContext = generateSingularityContext(indices, values, count)
            % Generate mask identifying singular interactions between rooftops.
            indexMask = rooftopMoM.core.impedanceMatrix.utils.generateSingularityMask(indices);

            % Compute unique rooftop size permutations and replicate over basis.
            [sizes, indicesPerm] = rooftopMoM.core.impedanceMatrix.utils.generateRooftopSizesUniquePerms(values);
            sizeMask = repmat(indicesPerm, 1, count);

            % Evaluate singular integral constants for each size case.
            type3Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type3(sizes(:, 1), sizes(:, 2));

            singularityContext = struct( ...
                'indexMask', indexMask, ...
                'sizeMask', sizeMask, ...
                'type3Values', type3Values ...
            );
        end

        function singularities = generateSingularities(singularityContext, idMask)
            % Define polarity groupings for same and opposite polarities.
            p1 = [1, 2, 3, 4];
            pCount = numel(p1);

            % Distribute singular integral values according to rooftop sizes.
            type3Values = singularityContext.type3Values(singularityContext.sizeMask);

            % Identify valid self-interactions from ID mask.
            validIndices = (idMask > 0);
            combinedMask = bsxfun(@and, singularityContext.indexMask, validIndices);

            % Allocate polarity-resolved singular contributions.
            singularities = combinedMask .* repmat(type3Values, [1, 1, pCount]);
        end
    end
end
