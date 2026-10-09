function [solution, converged, residualNorms] = iterate_(obj, apply, rightHandSide, initialGuess)
% iterate_: Runs the restarted GMRES, see rooftopMoM.solver.IterativeMethod.

    logger = rooftopMoM.logging.Logger.instance();
    isTrace = logger.isEnabled(rooftopMoM.logging.LogLevel.TRACE);
    count = numel(rightHandSide);
    cycleLength = min([obj.restart_, count, obj.maxIterationCount_]);

    reference = norm(obj.precondition_(rightHandSide));

    if isempty(initialGuess)
        solution = zeros(count, 1, 'like', 1i);
        residual = obj.precondition_(rightHandSide);
    else
        solution = initialGuess;
        residual = obj.precondition_(rightHandSide - apply(solution));
    end

    residualNorm = norm(residual);
    residualNorms = zeros(obj.maxIterationCount_ + 1, 1);
    residualNorms(1) = residualNorm / reference;
    converged = residualNorms(1) < obj.tolerance_;
    iteration = 0;

    basis = zeros(count, cycleLength + 1, 'like', 1i);
    hessenberg = zeros(cycleLength + 1, cycleLength, 'like', 1i);
    cosines = zeros(cycleLength, 1, 'like', 1i);
    sines = zeros(cycleLength, 1, 'like', 1i);

    while ~converged && iteration < obj.maxIterationCount_
        basis(:, 1) = residual / residualNorm;
        clear residual;

        % The right-hand side of the least squares problem, rotated along with the Hessenberg matrix.
        projection = [residualNorm; zeros(cycleLength, 1, 'like', 1i)];
        step = 0;

        while ~converged && step < cycleLength && iteration < obj.maxIterationCount_
            step = step + 1;
            iteration = iteration + 1;

            [vector, coefficients] = orthogonalize_(obj.precondition_(apply(basis(:, step))), basis, step);
            hessenberg(1:step, step) = coefficients(1:step);
            hessenberg(step + 1, step) = norm(vector);

            for k = 1:(step - 1)
                hessenberg([k, k + 1], step) = [conj(cosines(k)), conj(sines(k)); -sines(k), cosines(k)] ...
                    * hessenberg([k, k + 1], step);
            end

            [cosines(step), sines(step)] = generateRotation_(hessenberg(step, step), hessenberg(step + 1, step));
            subdiagonal = hessenberg(step + 1, step);
            hessenberg(step, step) = conj(cosines(step)) * hessenberg(step, step) + conj(sines(step)) * subdiagonal;
            hessenberg(step + 1, step) = 0;
            projection(step + 1) = -sines(step) * projection(step);
            projection(step) = conj(cosines(step)) * projection(step);

            residualNorms(iteration + 1) = abs(projection(step + 1)) / reference;
            converged = residualNorms(iteration + 1) < obj.tolerance_;

            if isTrace
                logger.trace('Iteration %d: relative residual %.3e.', iteration, residualNorms(iteration + 1));
            end

            % A zero subdiagonal means the basis spans the solution, and the residual is zero.
            if ~converged && subdiagonal ~= 0
                basis(:, step + 1) = vector / subdiagonal;
            end

            clear vector;
        end

        coefficients = zeros(cycleLength + 1, 1, 'like', 1i);
        coefficients(1:step) = hessenberg(1:step, 1:step) \ projection(1:step);
        solution = solution + basis * coefficients;

        % The next cycle starts from the true residual, as rounding lets the rotated one drift.
        if ~converged && iteration < obj.maxIterationCount_
            residual = obj.precondition_(rightHandSide - apply(solution));
            residualNorm = norm(residual);
        end
    end

    residualNorms = residualNorms(1:(iteration + 1));
end
