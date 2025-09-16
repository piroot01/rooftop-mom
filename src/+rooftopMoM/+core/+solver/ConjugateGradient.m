classdef ConjugateGradient < handle
% ConjugateGradient: MATLAB implementation of a (matrix-free) Conjugate Gradient solver
%   that multiplies vectors by a Green's-matrix operator implemented in
%   the spectral domain (FFT-based).

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
        function output = solve(obj, varargin)
            % solve: High-level entry point for the solver.
            %
            %   Optional name-value inputs:
            %       'initCurrent': Initial guess vector (complex, column).
            %       'excitation': Right-hand side vector (required).
            %       'tolerance': Relative residual tolerance (default 1e-6).
            %       'maxIterCount': Maximum number of iterations (default 1000).

            TRACE("ConjugateGradient.solve: Initialized solving procedure.");

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

            DEBUG("ConjugateGradient.solve: Parsed inputs: rooftop count = %d, tolerance = %g, maxIterCount = %d", ...
                count, p.Results.tolerance, p.Results.maxIterCount);

            output = obj.solveImpl_( ...
                p.Results.initCurrent, ...
                p.Results.excitation, ...
                p.Results.tolerance, ...
                p.Results.maxIterCount ...
            );
        end
    end

    methods (Access = protected)
        function output = solveImpl_(obj, initCurrent, excitation, tolerance, maxIterCount)
            % solveImpl_: Core Conjugate-Gradient loop.
            %
            %   Inputs:
            %       initCurrent: Initial guess column vector.
            %       excitation: Right-hand side column vector.
            %       tolerance: Relative tolerance.
            %       maxIterCount: Maximum iterations.
            %
            %   Returns:
            %       output: Solution in form of the ConjugateGradientOutput.

            TRACE("ConjugateGradient.solveImpl_: Starting CG solver core.");

            output = rooftopMoM.core.solver.ConjugateGradientOutput();

            current = initCurrent;
            residue = excitation - obj.applyGreensMatrix_(current);
            direction = residue;
            rTr = residue.' * residue;
            excitationNorm = norm(excitation);
            converged = false;

            DEBUG("ConjugateGradient.solveImpl_: Initial residual norm (unnormalized) = %g, excitationNorm = %g.", sqrt(rTr), excitationNorm);

            tic;

            for iteration = 1:maxIterCount
                DEBUG("ConjugateGradient.solveImpl_: Iteration %d: residual norm (prev) = %g", iteration, rTr);

                gDirection = obj.applyGreensMatrix_(direction);
                alpha = rTr / (direction.' * gDirection);
                current = current + alpha * direction;
                residue = residue - alpha * gDirection;
                rTrNew = residue.' * residue;
                residueNorm = sqrt(rTrNew) / excitationNorm;

                if residueNorm < tolerance
                    INFO("ConjugateGradient.solveImpl_: Converged in %d iterations (residue norm = %g).", iteration, residueNorm);
                    converged = true;
                    break;
                end

                beta = rTrNew / rTr;
                direction = residue + beta * direction;
                rTr = rTrNew;

                % Insert iteration info.
                output.insertIterationInfo( ...
                    'residueNorm', residueNorm ...
                );
            end

            solvingTime = toc;

            if ~converged
                INFO("ConjugateGradient.solveImpl_: Did not converge after %d iterations (final residue norm = %g).", maxIterCount, residueNorm, solvingTime);
            end

            INFO("ConjugateGradient.solveImpl_: Total time = %g s.", solvingTime);

            output.setPostRunInfo( ...
                struct( ...
                    'converged', converged, ...
                    'current', current, ...
                    'solvingTime', solvingTime, ...
                    'residueNorm', residueNorm ...
                ) ...
            );

            TRACE("ConjugateGradient.solveImpl_: Exiting CG solver core.");
        end

        function result = applyGreensMatrix_(obj, input)
            % applyGreensMatrix_: Compute y = A * input using FFT-based kernel.
            %   The operator depends on the rooftop basis content which can be
            %   XONLY, YONLY, or BOTH. The strategy in each case is:
            %   1) reshape/transpose input to 2D grid form matching the
            %      physical arrangement of rooftop basis functions;
            %   2) extend the grid to the padding size and compute fft2;
            %   3) multiply pointwise by the precomputed spectral kernel(s);
            %   4) inverse FFT (ifft2) to go back to spatial domain;
            %   5) truncate the padded result to the useful subregion and
            %      reshape back to a column vector if needed.

            sizeX = obj.context_.rooftops().countDir('x');
            sizeY = obj.context_.rooftops().countDir('y');

            switch obj.greensMatrix_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    % X-only basis functions: inputs are already in the
                    % expected vector layout (no explicit reshape here).
                    % Zero-extend and FFT.
                    inputFFT = fft2(input, obj.greensMatrix_.paddingX(), obj.greensMatrix_.paddingY());

                    % Spectral multiplication: G_xx * FFT(input).
                    outputFFT = obj.greensMatrix_.gxxFFT() .* inputFFT;

                    % Inverse FFT back to spatial domain.
                    output = ifft2(outputFFT);

                    % Truncate the useful part and return as a column (truncate returns an array with shape [sizeX, 1]).
                    result = rooftopMoM.core.solver.utils.truncate(output, sizeX, 1);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    % Y-only basis functions: the physical layout requires transposition before padding/FFT.
                    input = input.';

                    % Zero-extend and FFT (transpose handled above).
                    inputFFT = fft2(input, obj.greensMatrix_.paddingX(), obj.greensMatrix_.paddingY());

                    % Spectral multiplication.
                    outputFFT = obj.greensMatrix_.gyyFFT() .* inputFFT;

                    % Inverse FFT.
                    output = ifft2(outputFFT);

                    % Truncate the useful region and return vectorized.
                    outputTruncated = rooftopMoM.core.solver.utils.truncate(output, 1, sizeY);
                    result = outputTruncated(:);
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    % Both X and Y basis functions are present. The vector is split into two segments: 
                    % first the X-directed rooftop coefficients, then the Y-directed ones.
                    primitive = obj.context_.rooftops().primitive();
                    countX = primitive.countX();
                    countY = primitive.countY();

                    % Split vector into two parts corresponding to the two distinct grid samplings for X and Y roofops
                    inputX = input(1:sizeX);
                    inputY = input((sizeX + 1):end);

                    % Reshape into 2D arrays matching their grid layout.
                    inputXReshaped = reshape(inputX, countX - 1, countY);
                    inputYReshaped = reshape(inputY, countX, countY - 1);

                    % FFT each component after padding.
                    inputXFFT = fft2(inputXReshaped, obj.greensMatrix_.paddingX(), obj.greensMatrix_.paddingY());
                    inputYFFT = fft2(inputYReshaped, obj.greensMatrix_.paddingX(), obj.greensMatrix_.paddingY());

                    % Apply full 2x2 spectral Green's tensor.
                    outputXFFT = obj.greensMatrix_.gxxFFT() .* inputXFFT + obj.greensMatrix_.gxyFFT() .* inputYFFT;
                    outputYFFT = obj.greensMatrix_.gyxFFT() .* inputXFFT + obj.greensMatrix_.gyyFFT() .* inputYFFT;

                    % Inverse FFT back.
                    outputX = ifft2(outputXFFT);
                    outputY = ifft2(outputYFFT);

                    % Truncate to useful region (removing padding) and re-vectorize both parts to produce final column.
                    outputXTruncated = rooftopMoM.core.solver.utils.truncate(outputX, countX - 1, countY);
                    outputYTruncated = rooftopMoM.core.solver.utils.truncate(outputY, countX, countY - 1);

                    result = [outputXTruncated(:); outputYTruncated(:)];
                otherwise
                    error('Unknown rooftop content: %s', rooftopContent);
            end
        end
    end
end
