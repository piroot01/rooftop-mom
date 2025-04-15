classdef ConjugateGradientSolver < handle
    properties (Access = private)
        strategyContext_;

        solverContext_;

        magPartXFFT_;
        magPartYFFT_;
        elPartFFT_;
    end

    methods
        function obj = ConjugateGradientSolver(strategyContext)
            obj.strategyContext_ = strategyContext;
        end
    end

    methods
        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrent(obj, excitationVectorX, excitationVectorY, tolerance, maxInterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            [currentVector, converged, residueNorm, iterationCount] = obj.conjugateGradientSolver(excitationVector, tolerance, maxInterationCount);

            currentVectorX = reshape(currentVector(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);
        end

        function generateOperators(obj)
            magPartStrategy = mom.FastImpedanceMatrix.ComputeMagPart();
            elPartStrategy = mom.FastImpedanceMatrix.ComputeElPart();

            [magPartX, magPartY] = magPartStrategy.compute(obj.strategyContext_);
            elPart = elPartStrategy.compute(obj.strategyContext_);

            gridSize = double(obj.strategyContext_.rooftopBasisFunctions.getGridSize());
            obj.magPartXFFT_ = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartX, gridSize));
            obj.magPartYFFT_ = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartY, gridSize));
            obj.elPartFFT_ = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(elPart, gridSize));

            obj.solverContext_.basisCountX.x = gridSize(1) - 1;
            obj.solverContext_.basisCountX.y = gridSize(2);
            obj.solverContext_.basisCountY.x = gridSize(1);
            obj.solverContext_.basisCountY.y = gridSize(2) - 1;
            obj.solverContext_.offsetX = (gridSize(1) - 1) * gridSize(2);

            % does not matter which operator is used
            [M, N] = size(obj.magPartXFFT_);
            [u, v] = ndgrid(0:(M - 1), 0:(N - 1));
            obj.solverContext_.forwardShiftX = exp((1i * pi * u) / gridSize(1));
            obj.solverContext_.forwardShiftY = exp((1i * pi * v) / gridSize(2));
            obj.solverContext_.backwardShiftX = conj(obj.solverContext_.forwardShiftX);
            obj.solverContext_.backwardShiftY = conj(obj.solverContext_.forwardShiftY);

            obj.solverContext_.paddingX = 2 * gridSize(1);
            obj.solverContext_.paddingY = 2 * gridSize(2);
        end
    end

    methods (Access = private)
        function result = applyOperators(obj, currentVector)
            currentVectorX = reshape(currentVector(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);

            paddingX = obj.solverContext_.paddingX;
            paddingY = obj.solverContext_.paddingY;
            currentVectorXFFT = fft2(currentVectorX, paddingX, paddingY);
            currentVectorYFFT = fft2(currentVectorY, paddingX, paddingY);

            GxxTerm = obj.magPartXFFT_ + (1 - obj.solverContext_.forwardShiftX) .* (1 - obj.solverContext_.backwardShiftX) .* obj.elPartFFT_;
            GxyTerm = (1 - obj.solverContext_.forwardShiftX) .* (1 - obj.solverContext_.backwardShiftY) .* obj.elPartFFT_;
            GyxTerm = (1 - obj.solverContext_.forwardShiftY) .* (1 - obj.solverContext_.backwardShiftX) .* obj.elPartFFT_;
            GyyTerm = obj.magPartYFFT_ + (1 - obj.solverContext_.forwardShiftY) .* (1 - obj.solverContext_.backwardShiftY) .* obj.elPartFFT_;

            Gx = ifft2(GxxTerm .* currentVectorXFFT) + ifft2(GxyTerm .* currentVectorYFFT);
            Gy = ifft2(GyxTerm .* currentVectorXFFT) + ifft2(GyyTerm .* currentVectorYFFT);

            GxSection = Gx(1:obj.solverContext_.basisCountX.x, 1:obj.solverContext_.basisCountX.y);
            GySection = Gy(1:obj.solverContext_.basisCountY.x, 1:obj.solverContext_.basisCountY.y);

            result = [GxSection(:); GySection(:)];
        end

        function [currentVector, converged, residueNorm, iterationCount] = conjugateGradientSolver(obj, excitationVector, tolerance, maxInterationCount)
            currentVector = zeros(size(excitationVector));

            residue = excitationVector - obj.applyOperators(currentVector);
            direction = residue;
            rTr = residue.' * residue;
            excitationNorm = norm(excitationVector);
            converged = false;

            for iterationCount = 1:maxInterationCount
                GDirection = obj.applyOperators(direction);
                alpha = rTr / (direction.' * GDirection);
                currentVector = currentVector + alpha * direction;
                residue = residue - alpha * GDirection;
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
        end
    end
end
