classdef MagneticPartAssembler < rooftopMoM.core.impedanceMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_ rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner;
    end

    methods
        function obj = MagneticPartAssembler(assemblerChunkPlanner)
            obj.assemblerChunkPlanner_ = assemblerChunkPlanner;
        end
    end

    methods
        function Z = generateFull(obj, context)
            chunks = obj.assemblerChunkPlanner_.generateFull(context, true);

            Z = obj.generateSpecified(context, chunks);
        end

        function Z = generateSpecified(obj, context, chunks)
            wavenumber = context.wavenumber();
            matrixCount = numel(wavenumber);

            rooftops = context.rooftops();
            quadrature = context.quadrature();

            basisSize = rooftops.count();
            order = quadrature.order();
            blockSize = order^2;

            % Generate the quadrature, rooftop, rooftopDiv points (same for every matrix).
            [points, values] = context.rooftops().generateRooftops();
            [quadraturePoints, rooftopPoints] = rooftopMoM.core.impedanceMatrix.utils.convertRooftopPointsToArrays(points);

            slices = rooftopMoM.core.impedanceMatrix.utils.generateSlicesFromChunks(chunks, order);

            ids = rooftopMoM.basisFunctions.utils.generateIDMask(rooftops);
            idMask = rooftopMoM.core.impedanceMatrix.utils.generateIDMask(ids);
            jacobians = rooftopMoM.core.impedanceMatrix.utils.generate2DJacobians(rooftops.mesh(), ids);

            singularityContext = obj.generateSingularityContext_(rooftops.indices(), values, basisSize);

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

            % Pre-allocate the final impedance matrix (matrices).
            Z = complex(zeros(basisSize, basisSize, matrixCount));

            for n = 1:matrixCount
                Z(:, :, n) = obj.generateImpl_(computeContext, chunks, slices, wavenumber(n));
            end
        end
    end

    methods (Access = private)
        function Z = generateImpl_(obj, context, chunks, slices, wavenumber)
            basisSize = context.basisSize;
            chunkCount = numel(chunks);

            jacobians = context.jacobians;

            p1 = [1, 1, 2, 2];
            p2 = [1, 2, 1, 2];

            Z = complex(zeros(basisSize, basisSize));
            for n = 1:chunkCount
                % for debugging
                fprintf("Chunk ID=%d\n", n);

                chunk = chunks{n};
                slice = slices{n};

                rowSlice = slice.rowIndicesStart:slice.rowIndicesEnd;
                colSlice = slice.columnIndicesStart:slice.columnIndicesEnd;

                chunkSize = [numel(chunk.rowIndices), numel(chunk.columnIndices)];

                rowRooftopPoints = context.rooftopPoints(:, rowSlice, :);
                columnRooftopPoints = context.rooftopPoints(:, colSlice, :);

                singularityMask = context.singularityContext.indexMask(chunk.rowIndices, chunk.columnIndices, :);

                rowQuadraturePointsPolarity = context.quadraturePoints(:, rowSlice, p1);
                columnQuadraturePointsPolarity = context.quadraturePoints(:, colSlice, p2);

                distanceMatrix = rooftopMoM.core.impedanceMatrix.utils.generateDistanceMatrix(rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);
                integralKernel = rooftopMoM.core.impedanceMatrix.utils.generateScalarGreensFunction(wavenumber, distanceMatrix);

                idMask = context.idMask(chunk.rowIndices, chunk.columnIndices);
                sizeMask = context.singularityContext.sizeMask(chunk.rowIndices, chunk.columnIndices);

                chunkContext = struct( ...
                    'chunkSize', chunkSize, ...
                    'rowRooftopPoints', rowRooftopPoints, ...
                    'columnRooftopPoints', columnRooftopPoints, ...
                    'singularityMask', singularityMask, ...
                    'distanceMatrix', distanceMatrix, ...
                    'integralKernel', integralKernel, ...
                    'idMask', idMask, ...
                    'sizeMask', sizeMask ...
                );

                [matrixChunk, singularityMaskChunkNonZeroDot] = obj.generateChunk_(context, chunkContext, wavenumber);

                Z(chunk.rowIndices, chunk.columnIndices) = matrixChunk + obj.generateSingularities_(context.singularityContext, chunkContext, singularityMaskChunkNonZeroDot);
            end

            Z = 1i * wavenumber * rooftopMoM.core.Constants.Z0 / (4 * pi) .* Z .* (jacobians.' * jacobians);
        end

        function [chunk, singularityMaskChunkNonZeroDot] = generateChunk_(~, context, chunkContext, wavenumber)
            blockSize = context.blockSize;

            weights = repmat(context.weights, chunkContext.chunkSize);

            distanceMatrix = chunkContext.distanceMatrix;
            integralKernel = chunkContext.integralKernel;

            p1 = [1, 1, 2, 2];
            p2 = [1, 2, 1, 2];
            pCount = numel(p1);

            rowRooftopPointsPolarity = chunkContext.rowRooftopPoints(:, :, p1);
            columnRooftopPointsPolarity = chunkContext.columnRooftopPoints(:, :, p2);

            bigSingularityMask = repelem(chunkContext.singularityMask, blockSize, blockSize);

            distanceMatrixZero = (distanceMatrix < eps);

            rooftopPointsDot = pagemtimes(permute(rowRooftopPointsPolarity, [2, 1, 3]), columnRooftopPointsPolarity);
            rooftopPointsDotZero = (rooftopPointsDot < eps);

            bigSingularityMaskNonZeroDot = ~rooftopPointsDotZero & bigSingularityMask;
            singularityMaskChunkNonZeroDot = bigSingularityMaskNonZeroDot(1:blockSize:end, 1:blockSize:end, :);

            integralKernel(bigSingularityMaskNonZeroDot) = integralKernel(bigSingularityMaskNonZeroDot) - 1 ./ distanceMatrix(bigSingularityMaskNonZeroDot);

            integrands = rooftopPointsDot .* integralKernel;
            integrands(distanceMatrixZero) = rooftopPointsDot(distanceMatrixZero) * (-1i) * wavenumber;

            integrals = integrands .* weights;
            integralsReshaped = reshape(integrals, blockSize, chunkContext.chunkSize(1), blockSize, chunkContext.chunkSize(2), pCount);
            integralsSummed = sum(sum(integralsReshaped, 1), 3);
            chunkPolarity = reshape(integralsSummed, [chunkContext.chunkSize, pCount]);

            chunk = sum(chunkPolarity, 3);
        end

        function singularityContext = generateSingularityContext_(~, indices, values, count)
            indexMask = rooftopMoM.core.impedanceMatrix.utils.generateSingularityMask(indices);
            [sizes, indices] = rooftopMoM.core.impedanceMatrix.utils.generateRooftopSizesUniquePerms(values);

            sizeMask = repmat(indices, 1, count);
            type1Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type1(sizes(:, 1), sizes(:, 2));
            type2Values = rooftopMoM.core.impedanceMatrix.singularIntegrals.type2(sizes(:, 1), sizes(:, 2));

            singularityContext = struct( ...
                'indexMask', indexMask, ...
                'sizeMask', sizeMask, ...
                'type1Values', type1Values, ...
                'type2Values', type2Values ...
            );
        end

        function singularities = generateSingularities_(~, singularityContext, chunkContext, singularityMask)
            samePolarities = [1, 4];
            differentPolarities = [2, 3];
            pCount = numel(samePolarities) + numel(differentPolarities);

            sizeMask = chunkContext.sizeMask;
            type1Values = singularityContext.type1Values(sizeMask);
            type2Values = singularityContext.type2Values(sizeMask);

            validIndices = (chunkContext.idMask > 0);
            combinedMask = bsxfun(@and, singularityMask, validIndices);

            singularitiesPolarity = zeros([chunkContext.chunkSize, pCount]);

            singularitiesPolarity(:, :, samePolarities) = combinedMask(:, :, samePolarities) .* repmat(type1Values, [1, 1, 2]);
            singularitiesPolarity(:, :, differentPolarities) = combinedMask(:, :, differentPolarities) .* repmat(type2Values, [1, 1, 2]);

            singularities = sum(singularitiesPolarity, 3);
        end
    end
end
