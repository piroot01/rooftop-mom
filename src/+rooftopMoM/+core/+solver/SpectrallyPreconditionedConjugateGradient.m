classdef SpectrallyPreconditionedConjugateGradient < handle
    properties (Access = protected)
        matrixApplicator_;
    end

    methods
        function obj = SpectrallyPreconditionedConjugateGradient(matrixApplicator)
            obj.matrixApplicator_ = matrixApplicator;
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

            [inputSizeX, inputSizeY] = obj.matrixApplicator_.inputSize();

            p = inputParser();

            addParameter(p, 'initCurrent', complex(zeros(inputSizeX, inputSizeY)));
            addParameter(p, 'excitation', []);
            addParameter(p, 'tolerance', 1e-6);
            addParameter(p, 'maxIterCount', 1000);

            parse(p, varargin{:});

            if isempty(p.Results.excitation)
                error('The excitation is empty!');
            end

            DEBUG("SpectrallyPreconditionedConjugateGradient.solve: Parsed inputs: rooftop count = %d, tolerance = %g, maxIterCount = %d", ...
                inputSizeX * inputSizeY, p.Results.tolerance, p.Results.maxIterCount);

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
            residue = excitation - obj.matrixApplicator_.apply(current);
            preconditionedResidue = obj.matrixApplicator_.applyPreconditioner(residue);
            direction = preconditionedResidue;
            rTr = residue.' * preconditionedResidue;
            excitationNorm = norm(excitation);
            converged = false;

            DEBUG("SpectrallyPreconditionedConjugateGradient.solveImpl_: Initial residual norm (unnormalized) = %g, excitationNorm = %g.", sqrt(rTr), excitationNorm);

            t0 = tic;

            for iteration = 1:maxIterCount
                t1 = tic;

                gDirection = obj.matrixApplicator_.apply(direction);
                alpha = rTr / (direction.' * gDirection);
                current = current + alpha * direction;
                residue = residue - alpha * gDirection;
                preconditionedResidue = obj.matrixApplicator_.applyPreconditioner(residue);
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

                iterationTime = toc(t1);
                TRACE("ConjugateGradient.solveImpl_: Iteration %d: residual norm = %g, iteration time = %g s.", iteration, rTr, iterationTime);
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
    end
end
