classdef (Abstract) Product
% Product: Base of the matrix-free products with a matrix of a Problem, rooftopMoM.Operator and
%   rooftopMoM.Preconditioner. The matrix acts on the rooftop coefficients, the x-directed rooftops first.

    properties (Access = protected)
        problem_ rooftopMoM.Problem;

        % rooftopMoM.assembler.greensMatrix.Applicator doing the product.
        applicator_;
    end

    methods
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
end
