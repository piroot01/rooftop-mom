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
        function Z = generate(obj, context)
            chunks = obj.assemblerChunkPlanner_.generateFull(context, true);

            Z = obj.generateChunks(context, chunks);
        end

        function Z = generateChunks(obj, context, chunks)
            wavenumber = context.wavenumber();
            matrixCount = numel(wavenumber);
            basisSize = context.rooftops().count();
            order = context.quadrature().order();
            blockSize = order^2;

            % Generate the quadrature, rooftop, rooftopDiv points (same for every matrix).
            [points, values] = context.rooftops().generateRooftops();
            [quadraturePoints, rooftopPoints] = rooftopMoM.core.impedanceMatrix.utils.convertRooftopPointsToArrays(points);

            slices = rooftopMoM.core.utils.generateSlicesFromChunks(chunks, order);

            idMask = rooftopPoints.basisFunctions.utils.generateIDMask(context.rooftops());
            jacobians = rooftopMoM.core.utils.generate2DJacobians(rooftops.mesh(), IDs);

            singularityContext = obj.generateSingularityContext_(context.rooftops().indices(), values, basisSize);

            computeContext = struct( ...
                'basisSize', basisSize, ...
                'blockSize', blockSize, ...
                'chunks', chunks, ...
                'slices', slices, ...
                'quadraturePoints', quadraturePoints, ...
                'rooftopPoints', rooftopPoints, ...
                'jacobians', jacobians, ...
                'idMask', idMask, ...
                'singularityContext', singularityContext ...
            );

            % Pre-allocate the final impedance matrix (matrices).
            Z = complex(zeros(basisSize, basisSize, matrixCount));

            for n = 1:matrixCount
                Z(:, :, n) = obj.generateImpl_(computeContext, wavenumber(n));
            end
        end
    end

    methods (Access = private)
        function Z = generateImpl_(obj, context, wavenumber)
            basisSize = context.basisSize;
            blockSize = context.blockSize;
            chunkCount = numel(context.chunks);

            p1 = [1, 1, 2, 2];
            p2 = [1, 2, 1, 2];
            pCount = numel(p1);

            Z = complex(zeros(basisSize, basisSize));
            for n = 1:chunkCount
                % for debugging
                fprintf("Chunk ID=%d\n", n);

                chunk = context.chunks{n};
                slice = context.slices{n};

                rowSlice = slice.rowIndicesStart:slice.rowIndicesEnd;
                colSlice = slice.columnIndicesStart:slice.columnIndicesEnd;

                chunkSize = [numel(chunk.rowIndices), numel(chunk.columnIndices)];

                rowRooftopPoints = context.rooftopPoints(:, rowSlice, :);
                columnRooftopPoints = context.rooftopPoints(:, colSlice, :);

                singularityMask = context.singularityContext.indexMask(chunk.rowIndices, chunk.columnIndices, :);

                rowQuadraturePointsPolarity = quadraturePoints(:, rowSlice, p1);
                columnQuadraturePointsPolarity = quadraturePoints(:, colSlice, p2);

                integralKernel = obj.generateIntegralKernel(wavenumber, rowQuadraturePointsPolarity, columnQuadraturePointsPolarity)
            end
        end

        function singularityContext = generateSingularityContext_(~, indices, values, count)
            indexMask = rooftopMoM.core.impedanceMatrix.utils.generateSingularityMask(indices);
            [sizes, indices] = rooftopMoM.core.utils.generateRooftopSizesUniquePerms(values);

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
    end
end
