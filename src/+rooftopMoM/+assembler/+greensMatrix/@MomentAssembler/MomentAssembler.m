classdef MomentAssembler < rooftopMoM.assembler.greensMatrix.Assembler
% MomentAssembler: Assembles the Green's matrix generators from seven moments of the Green's function over every pair
%   of cells (the near-field note in docs/nearfield). The rooftops are linear on their two cells, so one integration
%   per cell offset serves all three generators (utils.accumulateMoments). The pairs whose cells are less than three
%   longer cell sides apart use the near quadrature order, the others the far one (utils.generateReducedRule). On the
%   touching pairs the quadrature integrates the remainder of the singularity extraction, and the moments of the
%   extracted part are added in closed form (singularIntegrals.generateExtractedMoments).

    properties (Constant, Access = private)
        % Gap between the cells of a near pair in both directions, in longer cell sides.
        nearGap_ = 3;
    end

    methods
        % generate: Returns the Generators of an AssemblerContext.
        generators = generate(obj, context);
    end

    methods (Static)
        function bytes = estimateWorkingMemory(rooftops, orders)
        % estimateWorkingMemory: Returns the bytes of the largest chunk the assembly of the rooftops with the
        %   quadrature orders [near, far] can use: its near or its far pass in one chunk. A larger chunk changes
        %   nothing.

            primitive = rooftops.primitive();
            counts = [primitive.countX(), primitive.countY()];
            sizes = [primitive.sizeX(), primitive.sizeY()];
            nearCounts = countNearOffsets_(counts, sizes, rooftopMoM.assembler.greensMatrix.MomentAssembler.nearGap_);
            nearRule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(orders(1), ...
                round(sizes / min(sizes)));
            farRule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(orders(2), [1, 1]);

            % A pass runs over every offset of its range, the far one over the near offsets too, and a chunk takes 24
            % bytes per kernel value (integrateCellPairs_).
            bytes = 24 * max(prod(nearCounts) * numel(nearRule.nodesX) * numel(nearRule.nodesY), ...
                prod(counts) * numel(farRule.nodesX) * numel(farRule.nodesY));
        end
    end
end
