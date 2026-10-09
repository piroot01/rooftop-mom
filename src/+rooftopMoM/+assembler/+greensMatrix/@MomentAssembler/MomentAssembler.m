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
end
