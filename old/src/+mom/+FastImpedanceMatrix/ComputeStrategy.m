classdef (Abstract) ComputeStrategy < handle
    properties (Access = private, Constant)
        bytesPerElement_ = double(8);
        bytesInMegabyte_ = double(1e6);
    end

    methods
        function mask = getSingularityMask_(~, strategyContext, direction1, direction2, ID)
            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);
            polaritiesCount = strategyContext.constants.polaritiesCount;

            basis1 = strategyContext.basisFunctions.getFacesInDirection(direction1);
            basisFunctionFaces1 = [basis1{ID}];

            basis2 = strategyContext.basisFunctions.getFacesInDirection(direction2);
            basisFunctionFaces2 = [basis2{ID}];

            mask = (reshape(basisFunctionFaces1(p1, :).', size(basisFunctionFaces1, 2), 1, polaritiesCount) == reshape(basisFunctionFaces2(p2, :).', 1, size(basisFunctionFaces2, 2), polaritiesCount));
        end

        function mask = getIDMask(~, strategyContext)
            basisFuncIDs = double(strategyContext.basisFunctions.getIDs());
            mask = (basisFuncIDs == basisFuncIDs.') .* (basisFuncIDs.' * ones(1, numel(basisFuncIDs)));
        end

        function jacobians = getJacobians(~, strategyContext)
            rooftopSizes = strategyContext.rooftopBasisFunctions.getRooftopSizes();
            jacobians2D = prod(rooftopSizes, 2);
            jacobians = jacobians2D(:) * jacobians2D(:).';
        end

        function mask = getSingularityMask(~, strategyContext, direction1, direction2)
            basisFunctionCount = strategyContext.rooftopBasisFunctions.getCount();
            mask = false(basisFunctionCount, 1, 4);

            if direction1 == 'x' && direction2 == 'x'
                % first plus rooftop x first plus rooftop
                mask(1, 1, 1) = true;

                % first minus rooftop x first plus rooftop
                mask(2, 1, 3) = true;

                % first minus rooftop x first minus rooftop
                mask(1, 1, 4) = true;
            elseif direction1 == 'y' && direction2 == 'y'
                gridSize = strategyContext.rooftopBasisFunctions.getGridSize();
                % first plus rooftop x first plus rooftop
                mask(1, 1, 1) = true;

                % first minus rooftop x first plus rooftop
                mask(gridSize(1) + 1, 1, 3) = true;

                % first minus rooftop x first minus rooftop
                mask(1, 1, 4) = true;
            end
        end

        function [chunks, slices] = getChunks(obj, strategyContext, basisFunctionCount ,memory_MB)
            elementCount = memory_MB * obj.bytesInMegabyte_ / obj.bytesPerElement_;
            quadratureOrder = strategyContext.quadrature.getOrder;
            elementCountPerBlockSide = quadratureOrder^2;
            elementCountSide = elementCount / elementCountPerBlockSide;

            if elementCountSide < elementCountPerBlockSide
                error('The given memory does not suffice!');
            end

            blockCountPerChunk = floor(double(elementCountSide) / elementCountPerBlockSide);
            chunkCount = floor(basisFunctionCount / blockCountPerChunk);

            if mod(basisFunctionCount, blockCountPerChunk) ~= 0
                chunkCount = chunkCount + 1;
            end

            chunks = cell(chunkCount, 1);
            slices = cell(chunkCount, 1);

            for chunkIndex = 1:chunkCount
                chunkStart = (chunkIndex - 1) * blockCountPerChunk + 1;
                chunkEnd = min(chunkIndex * blockCountPerChunk, basisFunctionCount);
                chunks{chunkIndex} = chunkStart:chunkEnd;

                slicesStart = (chunkIndex - 1) * (blockCountPerChunk * elementCountPerBlockSide) + 1;
                slicesEnd = min(chunkIndex * (elementCountPerBlockSide * blockCountPerChunk), basisFunctionCount * elementCountPerBlockSide);
                slices{chunkIndex} = [slicesStart, slicesEnd];
            end
        end
    end
end
