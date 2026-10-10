classdef MomentAssembler < rooftopMoM.assembler.greensMatrix.Assembler
% MomentAssembler: Assembles the Green's matrix generators from seven moments of the Green's function over every pair
%   of cells (the near-field note in docs/nearfield). The rooftops are linear on their two cells, so one integration
%   per cell offset serves all three generators (utils.accumulateMoments). The pairs whose cells are less than three
%   longer cell sides apart use the near quadrature order, the others the far one (utils.generateReducedRule). The far
%   pairs share the points of their rule between neighbouring offsets, so their kernel is evaluated once on a lattice
%   of these points (integrateLattice_). On the touching pairs the quadrature integrates the remainder of the
%   singularity extraction, and the moments of the extracted part are added in closed form
%   (singularIntegrals.generateExtractedMoments).

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
        %   quadrature orders [near, far] can use: its near pass in one chunk or its far pass in one tile. A larger
        %   chunk changes nothing.

            primitive = rooftops.primitive();
            counts = [primitive.countX(), primitive.countY()];
            sizes = [primitive.sizeX(), primitive.sizeY()];
            nearCounts = countNearOffsets_(counts, sizes, rooftopMoM.assembler.greensMatrix.MomentAssembler.nearGap_);
            nearRule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(orders(1), ...
                round(sizes / min(sizes)));

            % The near pass takes 24 bytes per kernel value of each of its offsets (integrateCellPairs_). The far pass
            % runs over the lattice blocks 0 to countX times 0 to countY, the near offsets included, at the peak bytes
            % of a block (integrateLattice_).
            bytes = max(24 * prod(nearCounts) * numel(nearRule.nodesX) * numel(nearRule.nodesY), ...
                estimateBlockBytes_(orders(2)) * prod(counts + 1));
        end
    end
end
