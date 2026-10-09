classdef (Abstract) Product < handle
% Product: Base of the matrix-free products with a matrix of a Problem, rooftopMoM.Operator and
%   rooftopMoM.Preconditioner. The matrix acts on the rooftop coefficients, the x-directed rooftops first. The memory of
%   its spectra is held in the budget of rooftopMoM.Resources while the object lives.

    properties (Access = protected)
        problem_ rooftopMoM.Problem;

        % rooftopMoM.assembler.greensMatrix.Applicator doing the product.
        applicator_;
    end

    properties (Access = private)
        % Bytes held in the budget, so that delete returns exactly them.
        heldBytes_ double = 0;
    end

    methods
        function delete(obj)
        % delete: Returns the memory of the spectra to the budget.

            rooftopMoM.resources.Budget.instance().release(obj.heldBytes_);
        end

        function output = apply(obj, input)
        % apply: Returns the product of the matrix with the Nx1 vector of rooftop coefficients input.

            count = obj.problem_.rooftopCount();

            if ~isnumeric(input) || ~isequal(size(input), [count, 1])
                error('rooftopMoM:product:sizeMismatch', 'The input must be a numeric %dx1 vector.', count);
            end

            output = obj.applicator_.apply(input);
        end

        function problem = problem(obj)
            problem = obj.problem_;
        end

        function applicator = applicator(obj)
        % applicator: Returns the core product. It skips the checks of apply, so the solvers iterate with it.

            applicator = obj.applicator_;
        end

        function bytes = memory(obj)
        % memory: Returns the bytes of the stored spectra.

            bytes = obj.applicator_.memory();
        end

        function bytes = workingMemory(obj)
        % workingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included.

            bytes = obj.applicator_.workingMemory();
        end
    end

    methods (Access = protected)
        function requireMemory_(~, problem, stored, working, description)
        % requireMemory_: Raises an error if the spectra and the working memory of the product do not fit the memory,
        %   together with the generators if the problem has not assembled them yet.
        %
        %   Inputs:
        %       problem: rooftopMoM.Problem.
        %       stored: Bytes of the spectra.
        %       working: Bytes alive at the peak of one product.
        %       description: Character vector naming the product.

            pending = 0;

            if problem.memory() == 0
                pending = rooftopMoM.resources.estimateGenerators(problem.rooftops());
            end

            toText = @rooftopMoM.resources.utils.formatBytes;
            rooftopMoM.resources.Budget.instance().require(pending + stored + working, sprintf(['%s (spectra %s, ' ...
                'working %s, generators to assemble %s)'], description, toText(stored), toText(working), toText(pending)));
        end

        function holdMemory_(obj)
        % holdMemory_: Holds the memory of the built spectra in the budget.

            obj.heldBytes_ = obj.applicator_.memory();
            rooftopMoM.resources.Budget.instance().hold(obj.heldBytes_);
        end
    end
end
