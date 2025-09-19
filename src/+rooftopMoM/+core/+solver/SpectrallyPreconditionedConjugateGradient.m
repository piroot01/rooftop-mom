classdef SpectrallyPreconditionedConjugateGradient < rooftopMoM.core.solver.ConjugateGradient
    methods
        function obj = SpectrallyPreconditionedConjugateGradient(context, greensMatrix)
            obj@rooftopMoM.core.solver.ConjugateGradient(context, greensMatrix);
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

            TRACE("SpectrallyPreconditionedConjugateGradient.solve: Initialized solving procedure.");

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

            DEBUG("SpectrallyPreconditionedConjugateGradient.solve: Parsed inputs: rooftop count = %d, tolerance = %g, maxIterCount = %d", ...
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
            % solveImpl_: Core Spectrally Preconditioned Conjugate-Gradient loop.
            %
            %   Inputs:
            %       initCurrent: Initial guess column vector.
            %       excitation: Right-hand side column vector.
            %       tolerance: Relative tolerance.
            %       maxIterCount: Maximum iterations.
            %
            %   Returns:
            %       output: Solution in form of the ConjugateGradientOutput.

            TRACE("SpectrallyPreconditionedConjugateGradient.solveImpl_: Starting CG solver core.");

            output = rooftopMoM.core.solver.ConjugateGradientOutput();

            current = initCurrent;
            residue = excitation - obj.applyGreensMatrix_(current);
            preconditionedResidue = obj.applyPreconditioner_(residue);
            direction = preconditionedResidue;
            rTr = residue.' * preconditionedResidue;
            excitationNorm = norm(excitation);
            converged = false;

            DEBUG("SpectrallyPreconditionedConjugateGradient.solveImpl_: Initial residual norm (unnormalized) = %g, excitationNorm = %g.", sqrt(rTr), excitationNorm);

            t0 = tic;

            for iteration = 1:maxIterCount
                TRACE("SpectrallyPreconditionedConjugateGradient.solveImpl_: Iteration %d: residual norm (prev) = %g", iteration, rTr);

                gDirection = obj.applyGreensMatrix_(direction);
                alpha = rTr / (direction.' * gDirection);
                current = current + alpha * direction;
                residue = residue - alpha * gDirection;
                preconditionedResidue = obj.applyPreconditioner_(residue);
                rTrNew = residue.' * preconditionedResidue;
                residueNorm = sqrt(rTrNew) / excitationNorm;

                if residueNorm < tolerance
                    DEBUG("SpectrallyPreconditionedConjugateGradient.solveImpl_: Converged in %d iterations (residue norm = %g).", iteration, residueNorm);
                    converged = true;
                    break;
                end

                beta = rTrNew / rTr;
                direction = preconditionedResidue + beta * direction;
                rTr = rTrNew;

                % Insert iteration info.
                output.insertIterationInfo( ...
                    'residueNorm', residueNorm ...
                );
            end

            solvingTime = toc(t0);

            if ~converged
                DEBUG("SpectrallyPreconditionedConjugateGradient.solveImpl_: Did not converge after %d iterations (final residue norm = %g).", maxIterCount, residueNorm, solvingTime);
            end

            DEBUG("SpectrallyPreconditionedConjugateGradient.solveImpl_: Total time = %g s.", solvingTime);

            output.setPostRunInfo( ...
                struct( ...
                    'converged', converged, ...
                    'current', current, ...
                    'solvingTime', solvingTime, ...
                    'residueNorm', residueNorm ...
                ) ...
            );

            TRACE("SpectrallyPreconditionedConjugateGradient.solveImpl_: Exiting CG solver core.");
        end

        function result = applyPreconditioner_(obj, input)
        % applyPreconditioner_: Apply spectral preconditioner to a rooftop-coefficient vector.
        %
        %   Description: Compute an approximation of y = P^{-1} * x (the action of the spectral
        %       preconditioner P^{-1} on a vector x) using elementwise operations in
        %       the Fourier (spectral) domain. The method:
        %       1) rearranges / reshapes the input into the 2-D grid that matches the
        %           rooftop basis layout (when necessary),
        %       2) zero-extends (pads) the grid to the precomputed FFT size,
        %       3) computes fft2(input_padded),
        %       4) performs pointwise (Hadamard) spectral division by the Green's
        %           spectral components (or uses the 2x2 block inverse for BOTH),
        %       5) returns to spatial domain with ifft2,
        %       6) truncates the padded region and re-vectorizes to produce the final column.

            sizeX = obj.context_.rooftops().countDir('x');
            sizeY = obj.context_.rooftops().countDir('y');

            switch obj.greensMatrix_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    % X-only basis functions: inputs are already in the
                    % expected vector layout (no explicit reshape here).
                    % Zero-extend and FFT.
                    inputFFT = fft2(input, obj.greensMatrix_.paddingX(), obj.greensMatrix_.paddingY());

                    % Do the hadamard division in the spectra.
                    outputFFT = inputFFT ./ obj.greensMatrix_.gxxFFT;

                    % Inverse FFT back to spatial domain.
                    output = ifft2(outputFFT);

                    % Truncate the useful part and return as a column (truncate returns an array with shape [sizeX, 1]).
                    result = output(1:sizeX, 1);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    % Y-only basis functions: the physical layout requires transposition before padding/FFT.
                    input = input.';

                    % Zero-extend and FFT (transpose handled above).
                    inputFFT = fft2(input, obj.greensMatrix_.paddingX(), obj.greensMatrix_.paddingY());

                    % Do the hadamard division in the spectra.
                    outputFFT = inputFFT ./ obj.greensMatrix_.gyyFFT;

                    % Inverse FFT.
                    output = ifft2(outputFFT);

                    % Truncate the useful region and return vectorized.
                    outputTruncated = output(1, 1:sizeY);

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

                    % Apply the block inversion of the Green's matrix, essentially aplying the P^(-1).
                    outputXFFT = (obj.greensMatrix_.gyyFFT() .* inputXFFT - obj.greensMatrix_.gxyFFT() .* inputYFFT) ./ obj.greensMatrix_.determinant();
                    outputYFFT = (obj.greensMatrix_.gxxFFT() .* inputYFFT - obj.greensMatrix_.gyxFFT() .* inputXFFT) ./ obj.greensMatrix_.determinant();

                    % Inverse FFT back.
                    outputX = ifft2(outputXFFT);
                    outputY = ifft2(outputYFFT);

                    % Truncate to useful region (removing padding) and re-vectorize both parts to produce final column.
                    outputXTruncated = outputX(1:(countX - 1), 1:countY);
                    outputYTruncated = outputY(1:countX, 1:(countY - 1));

                    result = [outputXTruncated(:); outputYTruncated(:)];
                otherwise
                    error('Unknown rooftop content: %s', rooftopContent);
            end
        end
    end
end
