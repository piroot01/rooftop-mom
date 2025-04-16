classdef ConjugateGradientSolver < handle
    properties (Access = private)
        strategyContext_;

        solverContext_;
    end

    methods
        function obj = ConjugateGradientSolver(strategyContext)
            obj.strategyContext_ = strategyContext;
        end
    end

    methods
        function generateOperators(obj)
            magPartStrategy = mom.FastImpedanceMatrix.ComputeMagPart();
            elPartStrategy = mom.FastImpedanceMatrix.ComputeElPart();

            [magPartX, magPartY] = magPartStrategy.compute(obj.strategyContext_);
            elPart = elPartStrategy.compute(obj.strategyContext_);

            gridSize = double(obj.strategyContext_.rooftopBasisFunctions.getGridSize());
            magPartXFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartX, gridSize));
            magPartYFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartY, gridSize));
            elPartFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(elPart, gridSize));

            obj.solverContext_.basisCountX.x = gridSize(1) - 1;
            obj.solverContext_.basisCountX.y = gridSize(2);
            obj.solverContext_.basisCountY.x = gridSize(1);
            obj.solverContext_.basisCountY.y = gridSize(2) - 1;
            obj.solverContext_.offsetX = (gridSize(1) - 1) * gridSize(2);

            % does not matter which operator is used
            [M, N] = size(magPartXFFT);
            [u, v] = ndgrid(0:(M - 1), 0:(N - 1));
            forwardShiftX = exp((1i * pi * u) / gridSize(1));
            forwardShiftY = exp((1i * pi * v) / gridSize(2));
            backwardShiftX = conj(forwardShiftX);
            backwardShiftY = conj(forwardShiftY);

            obj.solverContext_.paddingX = 2 * gridSize(1);
            obj.solverContext_.paddingY = 2 * gridSize(2);

            obj.solverContext_.Gxx = magPartXFFT + (1 - forwardShiftX) .* (1 - backwardShiftX) .* elPartFFT;
            obj.solverContext_.Gxy = (1 - forwardShiftX) .* (1 - backwardShiftY) .* elPartFFT;
            obj.solverContext_.Gyx = (1 - forwardShiftY) .* (1 - backwardShiftX) .* elPartFFT;
            obj.solverContext_.Gyy = magPartYFFT + (1 - forwardShiftY) .* (1 - backwardShiftY) .* elPartFFT;
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrent(obj, excitationVectorX, excitationVectorY, tolerance, maxInterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            [currentVector, converged, residueNorm, iterationCount] = obj.conjugateGradientSolver(excitationVector, tolerance, maxInterationCount);
            %[currentVector, converged, residueNorm, iterationCount] = bicgstab(@(excitationVector) obj.applyOperators(excitationVector), excitationVector, tolerance, maxInterationCount);
            %[currentVector, converged, residueNorm, iterationCount] = cgs(@(excitationVector) obj.applyOperators(excitationVector), excitationVector, tolerance, maxInterationCount);

            currentVectorX = reshape(currentVector(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);
        end

        function impedanceMatrix = computeImpedanceMatrix1(obj)
            basisFunctionCount = obj.solverContext_.basisCountX.x * obj.solverContext_.basisCountX.y + obj.solverContext_.basisCountY.x * obj.solverContext_.basisCountY.y;
            impedanceMatrix = zeros(basisFunctionCount);

            for column = 1:basisFunctionCount
                canonicalVector = zeros(basisFunctionCount, 1);
                canonicalVector(column) = 1;

                impedanceMatrix(:, column) = obj.applyOperatorsVector(canonicalVector);
            end
        end

        function impedanceMatrix = computeImpedanceMatrix2(obj)
            basisFunctionCount = obj.solverContext_.basisCountX.x * obj.solverContext_.basisCountX.y + obj.solverContext_.basisCountY.x * obj.solverContext_.basisCountY.y;

            I = eye(basisFunctionCount);

            impedanceMatrix = obj.applyOperatorsMatrix(I);
        end
    end

    methods (Access = private)
        function result = applyOperatorsVector(obj, currentVector)
            currentVectorX = reshape(currentVector(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);

            paddingX = obj.solverContext_.paddingX;
            paddingY = obj.solverContext_.paddingY;

            currentVectorXFFT = fft2(currentVectorX, paddingX, paddingY);
            currentVectorYFFT = fft2(currentVectorY, paddingX, paddingY);

            Gx = ifft2(obj.solverContext_.Gxx .* currentVectorXFFT + obj.solverContext_.Gxy .* currentVectorYFFT);
            Gy = ifft2(obj.solverContext_.Gyx .* currentVectorXFFT + obj.solverContext_.Gyy .* currentVectorYFFT);

            GxSection = Gx(1:obj.solverContext_.basisCountX.x, 1:obj.solverContext_.basisCountX.y);
            GySection = Gy(1:obj.solverContext_.basisCountY.x, 1:obj.solverContext_.basisCountY.y);

            result = [GxSection(:); GySection(:)];
        end

        function result = applyOperatorsMatrix(obj, currentMatrix)
            bcx = obj.solverContext_.basisCountX;
            bcy = obj.solverContext_.basisCountY;
            padX = obj.solverContext_.paddingX;
            padY = obj.solverContext_.paddingY;

            N = size(currentMatrix, 2);

            currentX = reshape(currentMatrix(1:bcx.x * bcx.y, :), bcx.x, bcx.y, N);
            currentY = reshape(currentMatrix(bcx.x * bcx.y+1:end, :), bcy.x, bcy.y, N);

            currents = zeros(padX, padY, 2, N);
            currents(1:bcx.x, 1:bcx.y, 1, :) = currentX;
            currents(1:bcy.x, 1:bcy.y, 2, :) = currentY;

            currentsFFT = fft2(currents);

            GxFFT = obj.solverContext_.Gxx .* currentsFFT(:, :, 1, :) + obj.solverContext_.Gxy .* currentsFFT(:, :, 2, :);
            GyFFT = obj.solverContext_.Gyx .* currentsFFT(:, :, 1, :) + obj.solverContext_.Gyy .* currentsFFT(:, :, 2, :);

            Gx = ifft2(GxFFT);
            Gy = ifft2(GyFFT);

            GxSection = Gx(1:bcx.x, 1:bcx.y, 1, :);
            GySection = Gy(1:bcy.x, 1:bcy.y, 1, :);

            result = [reshape(GxSection, bcx.x * bcx.y, N); reshape(GySection, bcy.x * bcy.y, N)];
        end

        function [currentVector, converged, residueNorm, iterationCount] = conjugateGradientSolver(obj, excitationVector, tolerance, maxInterationCount)
            currentVector = randn(size(excitationVector));

            residue = excitationVector - obj.applyOperatorsVector(currentVector);
            direction = residue;
            rTr = residue.' * residue;
            excitationNorm = norm(excitationVector);
            converged = false;

            for iterationCount = 1:maxInterationCount
                GDirection = obj.applyOperatorsVector(direction);
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
