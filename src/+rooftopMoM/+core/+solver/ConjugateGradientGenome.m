classdef ConjugateGradientGenome < rooftopMoM.core.solver.ConjugateGradient
    methods
        function obj = ConjugateGradientGenome(context, greensMatrix)
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
        %       'genome': Masks the unused rooftops.

            TRACE("ConjugateGradientGenome.solve: Initialized solving procedure.");

            count = obj.context_.rooftops().count();

            p = inputParser();

            addParameter(p, 'initCurrent', []);
            addParameter(p, 'excitation', []);
            addParameter(p, 'tolerance', 1e-6);
            addParameter(p, 'maxIterCount', 1000);
            addParameter(p, 'genome', true(count, 1));

            parse(p, varargin{:});

            if isempty(p.Results.excitation)
                error('The excitation is empty!');
            end

            if isempty(p.Results.initCurrent)
                initCurrent = complex(zeros(sum(p.Results.genome), 1));
            else
                initCurrent = p.Results.initCurrent;
            end

            DEBUG("ConjugateGradientGenome.solve: Parsed inputs: rooftop count = %d, tolerance = %g, maxIterCount = %d", ...
                count, p.Results.tolerance, p.Results.maxIterCount);

            output = obj.solveImpl_( ...
                initCurrent, ...
                p.Results.excitation, ...
                p.Results.tolerance, ...
                p.Results.maxIterCount, ...
                p.Results.genome ...
            );
        end
    end

    methods (Access = protected)
        function output = solveImpl_(obj, initCurrent, excitation, tolerance, maxIterCount, genome)
        % solveImpl_: Core Conjugate-Gradient loop.
        %
        %   Inputs:
        %       initCurrent: Initial guess column vector.
        %       excitation: Right-hand side column vector.
        %       tolerance: Relative tolerance.
        %       maxIterCount: Maximum iterations.
        %       genome: Logical vector selecting the active rooftop basis functions indices.
        %
        %   Returns:
        %       output: Solution in form of the ConjugateGradientOutput.

            TRACE("ConjugateGradientGenome.solveImpl_: Starting CG solver core.");

            output = rooftopMoM.core.solver.ConjugateGradientOutput();

            current = initCurrent;
            residue = excitation - obj.applyGreensMatrix_(current, genome);
            direction = residue;
            rTr = residue.' * residue;
            excitationNorm = norm(excitation);
            converged = false;

            DEBUG("ConjugateGradientGenome.solveImpl_: Initial residual norm (unnormalized) = %g, excitationNorm = %g.", sqrt(rTr), excitationNorm);

            tic;

            for iteration = 1:maxIterCount
                DEBUG("ConjugateGradientGenome.solveImpl_: Iteration %d: residual norm (prev) = %g", iteration, rTr);

                gDirection = obj.applyGreensMatrix_(direction, genome);
                alpha = rTr / (direction.' * gDirection);
                current = current + alpha * direction;
                residue = residue - alpha * gDirection;
                rTrNew = residue.' * residue;
                residueNorm = sqrt(rTrNew) / excitationNorm;

                if residueNorm < tolerance
                    INFO("ConjugateGradientGenome.solveImpl_: Converged in %d iterations (residue norm = %g).", iteration, residueNorm);
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
                INFO("ConjugateGradientGenome.solveImpl_: Did not converge after %d iterations (final residue norm = %g).", maxIterCount, residueNorm, solvingTime);
            end

            INFO("ConjugateGradientGenome.solveImpl_: Total time = %g s.", solvingTime);

            output.setPostRunInfo( ...
                struct( ...
                    'converged', converged, ...
                    'current', current, ...
                    'solvingTime', solvingTime, ...
                    'residueNorm', residueNorm ...
                ) ...
            );

            TRACE("ConjugateGradientGenome.solveImpl_: Exiting CG solver core.");
        end

        function result = applyGreensMatrix_(obj, input, genome)
            extendedInput = zeros(numel(genome), 1, 'like', input);
            extendedInput(genome) = input;

            output = applyGreensMatrix_@rooftopMoM.core.solver.ConjugateGradient(obj, extendedInput);

            result = output(genome);
        end
    end
end
