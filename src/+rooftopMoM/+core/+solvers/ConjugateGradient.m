classdef ConjugateGradient < handle
    properties (Access = protected)
        context_;

        greensMatrix_ rooftopMoM.core.greensMatrix.MatrixData;
    end

    methods
        function obj = ConjugateGradient(context, greensMatrix)
            obj.context_ = context;
            obj.greensMatrix_ = greensMatrix;
        end
    end

    methods
        function result = solve(obj, varargin)
            count = obj.context_.rooftops().count();

            p = inputParser();

            addParameter(p, 'initCurrent', complex(zeros(count, 1)));
            addParameter(p, 'excitation', []);
            addParameter(p, 'tolerance', 1e-6);
            addParameter(p, 'maxIterCount', 1000);

            parse(p, varargin{:});

            if isempty(p.Results.excitation)
                error('The excitation is empty!');
            end

            result = obj.solveImpl_( ...
                p.Results.initCurrent, ...
                p.Results.excitation, ...
                p.Results.tolerance, ...
                p.Results.maxIterCount ...
            );
        end
    end

    methods (Access = protected)
        function result = solveImpl_(obj, initCurrent, excitation, tolerance, maxIterCount)
            current = initCurrent;
            residue = excitation - obj.applyGreensMatrix_(current);
            direction = residue;
            rTr = residue.' * residue;
            excitationNorm = norm(excitation);
            converged = false;

            for iteration = 1:maxIterCount
                gDirection = obj.applyGreensMatrix_(direction);
                alpha = rTr / (direction.' * gDirection);
                current = current + alpha * direction;
                residue = residue - alpha * gDirection;
                rTrNew = residue.' * residue;
                residueNorm = sqrt(rTrNew) / excitationNorm;

                if residueNorm < tolerance
                    converged = true;
                    break;
                end

                beta = rTrNew / rTr;
                direction = residue + beta * direction;
                rTr = rTrNew;
            end

            result = current;
        end

        function result = applyGreensMatrix_(obj, input)
            countX = obj.context_.rooftops().countDir('x');
            countY = obj.context_.rooftops().countDir('y');

            switch obj.greensMatrix_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    % First we reshape the input based on the antenna shape.
                    input = input.';

                    % Then we apply the extension operator and perform the FFT at one go.
                    inputFFT = fft2(input, obj.greensMatrix_.paddingX(), obj.greensMatrix_.paddingY());

                    % Apply the green's matrix in its spectral form.
                    outputFFT = obj.greensMatrix_.gyyFFT() .* inputFFT;

                    % Go back from the reciprocal coordinates.
                    output = ifft2(outputFFT);

                    % Use the truncation operator for the extraction of the useful part of the output.
                    outputTruncated = rooftopMoM.core.solvers.utils.truncate(output, 1, countY);

                    % And finally vectorize the outputTruncated
                    result = outputTruncated(:);
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                otherwise
                    error('Unknown rooftop content: %s', rooftopContent);
            end
        end
    end
end
