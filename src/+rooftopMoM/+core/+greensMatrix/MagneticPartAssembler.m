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

            singularityContext = rooftopMoM.core.greensMatrix.MagneticPartAssembler.generateSingularityContext(values, basisSize);

            weights = quadrature.weights4D();
            G = [];
        end
    end

    methods (Static)
        function singularityContext = generateSingularityContext(values, count)
            indices = 1:count

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
