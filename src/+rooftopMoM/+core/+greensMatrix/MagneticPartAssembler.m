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
            chunksX = obj.assemblerChunkPlanner_.generateFullDir(context, 'x');
            chunksY = obj.assemblerChunkPlanner_.generateFullDir(context, 'y');

            [Gxx, Gyy] = obj.generateSpecified(context, chunksX, chunksY);
        end

        function [Gxx, Gyy] = generateSpecified(obj, context, chunksX, chunksY)
            Gxx = obj.generateSpecifiedDir(context, chunksX, 'x');
            Gyy = obj.generateSpecifiedDir(context, chunksY, 'y');
        end

        function G = generateFullDir(obj, context, direction)
            % Produce all chunks for full matrix assembly
            chunks = obj.assemblerChunkPlanner_.generateFullDir(context, direction);

            G = obj.generateImpl_(context, chunks, direction);
        end

        function G = generateSpecifiedDir(obj, context, chunks, direction)
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

            basisSize = rooftops.countDir(direction);
            order = quadrature.order();
            blockSize = order^2;

            [rooftopPoints, rooftopValues] = rooftops.generateRooftopsDir(direction);

            %{
            % Cache the data
            if rooftopMoM.Cache.getStatus()
                quadraturePointsKey = strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_QUADRATUREPOINTS, upper(direction));
                rooftopPointsKey = strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPPOINTS, upper(direction));
                rooftopValuesKey = strcat(rooftopMoM.CacheKeys.CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPVALUES, upper(direction));

                rooftopMoM.Cache.autoSet(quadraturePointsKey, quadraturePoints);
                rooftopMoM.Cache.autoSet(rooftopPointsKey, rooftopPoints);
                rooftopMoM.Cache.autoSet(rooftopValuesKey, values);
            end
            %}

            jacobians = rooftopMoM.core.greensMatrix.utils.generate2DJacobians(rooftops.primitive(), basisSize);

            singularityContext = rooftopMoM.core.greensMatrix.MagneticPartAssembler.generateSingularityContext(rooftops.indicesDir(direction), rooftopValues);

            weights = quadrature.weights4D();

            computeContext = struct( ...
                'direction', direction, ...
                'primitive', rooftops.primitive(), ...
                'basisSize', basisSize, ...
                'blockSize', blockSize, ...
                'rooftopPoints', rooftopPoints, ...
                'jacobians', jacobians, ...
                'singularityContext', singularityContext, ...
                'weights', weights ...
            );

            % TODO: Include the multi-wavenumber output
            G = obj.generateImplCore_(computeContext, chunks, wavenumber);

            % Reshape according to the primitive shape.
            [cX, cY] = rooftopMoM.core.greensMatrix.utils.generateGreensMatrixShape(rooftops.primitive(), direction);

            % NOTE: the cY is aligned with columns and the cX with the rows
            G = reshape(G, cX, cY);
        end

        function G = generateImplCore_(~, context, chunks, wavenumber)
            basisSize = context.basisSize;
            jacobians = context.jacobians;
            % Initialize greens matrix for this frequency.
            G = complex(zeros(basisSize, 1));

            % Process each chunk sequentially.
            for n = 1:numel(chunks)
                DEBUG("MagneticPartAssembler.generateImplCore_: Assembling chunk %d/%d.", n, numel(chunks));

                chunk = chunks{n};
                chunkRowSize = numel(chunk.rowIndices);
                chunkColumnSize = numel(chunk.columnIndices);

                rowRooftopPointsPositive = repmat(context.rooftopPoints.rooftopPointsPositive, 1, chunkRowSize);
                rowRooftopPointsNegative = repmat(context.rooftopPoints.rooftopPointsNegative, 1, chunkRowSize);
                columnRooftopPointsPositive = repmat(context.rooftopPoints.rooftopPointsPositive, 1, chunkColumnSize);
                columnRooftopPointsNegative = repmat(context.rooftopPoints.rooftopPointsNegative, 1, chunkColumnSize);

                rowRooftopPoints = cat(3, rowRooftopPointsPositive, rowRooftopPointsNegative);
                columnRooftopPoints = cat(3, columnRooftopPointsPositive, columnRooftopPointsNegative);

                % Row quadrature points in-place generation.
                rowShifts = rooftopMoM.basisFunctions.utils.generateRectGridShiftsSliceDirection( ...
                    context.primitive.countX, context.primitive.sizeX, context.primitive.sizeY, ...
                    chunk.rowIndices, ...
                    context.direction ...
                );

                rowQuadraturePointsPositive = rooftopMoM.core.greensMatrix.utils.generateQuadraturePointsSlice(context.rooftopPoints.quadraturePointsPositive, rowShifts);
                rowQuadraturePointsNegative = rooftopMoM.core.greensMatrix.utils.generateQuadraturePointsSlice(context.rooftopPoints.quadraturePointsNegative, rowShifts);

                % Column quadrature points in-place generation.
                columnShifts = rooftopMoM.basisFunctions.utils.generateRectGridShiftsSliceDirection( ...
                    context.primitive.countX, context.primitive.sizeX, context.primitive.sizeY, ...
                    chunk.columnIndices, ...
                    context.direction ...
                );

                columnQuadraturePointsPositive = rooftopMoM.core.greensMatrix.utils.generateQuadraturePointsSlice(context.rooftopPoints.quadraturePointsPositive, columnShifts);
                columnQuadraturePointsNegative = rooftopMoM.core.greensMatrix.utils.generateQuadraturePointsSlice(context.rooftopPoints.quadraturePointsNegative, columnShifts);

                rowQuadraturePointsPolarity = cat(3, rowQuadraturePointsPositive, rowQuadraturePointsPositive, rowQuadraturePointsNegative, rowQuadraturePointsNegative);
                columnQuadraturePointsPolarity = cat(3, columnQuadraturePointsPositive, columnQuadraturePointsNegative, columnQuadraturePointsPositive, columnQuadraturePointsNegative);

                distanceMatrix = rooftopMoM.core.impedanceMatrix.utils.generateDistanceMatrix(rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);
                integralKernel = rooftopMoM.core.impedanceMatrix.utils.generateScalarGreensFunction(wavenumber, distanceMatrix);

                % Obtain mask of singular interactions.
                singularityMask = context.singularityContext.indexMask(chunk.rowIndices, chunk.columnIndices, :);

                % Bundle chunk-specific data for integration.
                chunkContext = struct( ...
                    'chunkSize', [chunkRowSize, chunkColumnSize], ...
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

            % Evaluate the singular integral constants for the given size.
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
