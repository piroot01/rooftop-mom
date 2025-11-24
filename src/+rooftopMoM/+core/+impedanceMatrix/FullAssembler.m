classdef FullAssembler < rooftopMoM.core.impedanceMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_ rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner;
    end

    methods
        function obj = FullAssembler(assemblerChunkPlanner)
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
            rooftopDivPoints = rooftopMoM.core.impedanceMatrix.utils.stackPositiveAndNegativeFields(values, 'rooftopDivPositive', 'rooftopDivNegative');

            % Cache points.
            rooftopMoM.Cache.autoSet(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS, rooftopPoints);

            slices = rooftopMoM.core.impedanceMatrix.utils.generateSlicesFromChunks(chunks, order);

            ids = rooftopMoM.basisFunctions.utils.generateIDMask(rooftops);
            idMask = rooftopMoM.core.impedanceMatrix.utils.generateIDMask(ids);
            jacobians = rooftopMoM.core.impedanceMatrix.utils.generate2DJacobians(rooftops.mesh(), ids);

            singularityContext = rooftopMoM.core.impedanceMatrix.FullAssembler.generateSingularityContext(rooftops.indices(), values, basisSize);
            electricSingularities = rooftopMoM.core.impedanceMatrix.ElectricPartAssembler.generateSingularities(singularityContext, idMask);

            weights = quadrature.weights4D();

            computeContext = struct( ...
                'basisSize', basisSize, ...
                'blockSize', blockSize, ...
                'quadraturePoints', quadraturePoints, ...
                'rooftopPoints', rooftopPoints, ...
                'rooftopDivPoints', rooftopDivPoints, ...
                'jacobians', jacobians, ...
                'idMask', idMask, ...
                'singularityContext', singularityContext, ...
                'electricSingularities', electricSingularities, ...
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
                DEBUG("FullAssembler.generateImpl_: Assembling chunk %d/%d.", n, numel(chunks));

                chunk = chunks{n};
                slice = slices{n};

                % Determine row and column index ranges for current chunk.
                rowSlice = slice.rowIndicesStart:slice.rowIndicesEnd;
                colSlice = slice.columnIndicesStart:slice.columnIndicesEnd;

                rowRooftopPoints = context.rooftopPoints(:, rowSlice, :);
                columnRooftopPoints = context.rooftopPoints(:, colSlice, :);

                rowRooftopDivPoints = context.rooftopDivPoints(:, chunk.rowIndices, :);
                columnRooftopDivPoints = context.rooftopDivPoints(:, chunk.columnIndices, :);

                rowQuadraturePointsPolarity = context.quadraturePoints(:, rowSlice, p1);
                columnQuadraturePointsPolarity = context.quadraturePoints(:, colSlice, p2);

                distanceMatrix = rooftopMoM.core.impedanceMatrix.utils.generateDistanceMatrix(rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);
                integralKernel = rooftopMoM.core.impedanceMatrix.utils.generateScalarGreensFunction(wavenumber, distanceMatrix);

                % Extract ID and size masks for singularity correction.
                idMask = context.idMask(chunk.rowIndices, chunk.columnIndices);
                sizeMask = context.singularityContext.sizeMask(chunk.rowIndices, chunk.columnIndices);

                % Obtain mask and values (eldctric part) of singular interactions.
                singularityMask = context.singularityContext.indexMask(chunk.rowIndices, chunk.columnIndices, :);
                electricSingularities = context.electricSingularities(chunk.rowIndices, chunk.columnIndices, :);

                % Bundle chunk-specific data for integration.
                chunkContext = struct( ...
                    'chunkSize', [numel(chunk.rowIndices), numel(chunk.columnIndices)], ...
                    'rowRooftopPoints', rowRooftopPoints, ...
                    'columnRooftopPoints', columnRooftopPoints, ...
                    'rowRooftopDivPoints', rowRooftopDivPoints, ...
                    'columnRooftopDivPoints', columnRooftopDivPoints, ...
                    'singularityMask', singularityMask, ...
                    'electricSingularities', electricSingularities, ...
                    'distanceMatrix', distanceMatrix, ...
                    'integralKernel', integralKernel, ...
                    'idMask', idMask, ...
                    'sizeMask', sizeMask ...
                );

                % Compute non-singular and singular contributions for the magnetic part.
                [magneticMatrixChunk, singularityMaskChunkNonZeroDot] = rooftopMoM.core.impedanceMatrix.MagneticPartAssembler.generateChunk(context, chunkContext, wavenumber);
                magneticSingularCorrection = rooftopMoM.core.impedanceMatrix.MagneticPartAssembler.generateSingularities(context.singularityContext, chunkContext, singularityMaskChunkNonZeroDot);

                % Apply the singular corrections.
                magneticMatrixChunk = magneticMatrixChunk + magneticSingularCorrection;

                % Compute constributions for the electric part.
                electricMatrixChunk = rooftopMoM.core.impedanceMatrix.ElectricPartAssembler.generateChunk(context, chunkContext, wavenumber);

                % Scale the electric part by the factor 1/k^2.
                electricMatrixChunk = electricMatrixChunk .* (1 / wavenumber^2);

                % Accumulate into the global Z matrix.
                Z(chunk.rowIndices, chunk.columnIndices) = magneticMatrixChunk - electricMatrixChunk;
            end

            % Apply final scaling by free-space impedance and Jacobian factors.
            Z = 1i * wavenumber * rooftopMoM.core.Constants.Z0 / (4 * pi) .* Z .* (jacobians.' * jacobians);
        end
    end

    methods (Static)
        function singularityContext = generateSingularityContext(indices, values, count)
            % Generate mask identifying singular interactions between rooftops.
            indexMask = rooftopMoM.core.impedanceMatrix.utils.generateSingularityMask(indices);

            % Compute unique rooftop size permutations and replicate over basis.
            [sizes, indicesPerm] = rooftopMoM.core.impedanceMatrix.utils.generateRooftopSizesUniquePerms(values);
            sizeMask = repmat(indicesPerm, 1, count);

            % Evaluate singular integral constants for each size case.
            type1Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type1(sizes(:, 1), sizes(:, 2));
            type2Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type2(sizes(:, 1), sizes(:, 2));
            type3Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type3(sizes(:, 1), sizes(:, 2));

            singularityContext = struct( ...
                'indexMask', indexMask, ...
                'sizeMask', sizeMask, ...
                'type1Values', type1Values, ...
                'type2Values', type2Values, ...
                'type3Values', type3Values ...
            );
        end
    end
end
