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

            obj.solverContext_.basisCountX.x = gridSize(1) - 1;
            obj.solverContext_.basisCountX.y = gridSize(2);
            obj.solverContext_.basisCountY.x = gridSize(1);
            obj.solverContext_.basisCountY.y = gridSize(2) - 1;
            obj.solverContext_.offsetX = (gridSize(1) - 1) * gridSize(2);

            M = gridSize(1) * 2;
            N = gridSize(2) * 2;

            obj.solverContext_.paddingX = M;
            obj.solverContext_.paddingY = N;

            [u, v] = ndgrid(0:(M - 1), 0:(N - 1));
            forwardShiftX = exp((1i * pi * u) / gridSize(1));
            forwardShiftY = exp((1i * pi * v) / gridSize(2));
            backwardShiftX = conj(forwardShiftX);
            backwardShiftY = conj(forwardShiftY);

            if ~all(isnan(magPartX), 'all') && all(isnan(magPartY), 'all')
                obj.solverContext_.rooftopContent = 'onlyX';

                magPartXFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartX, gridSize));
                elPartFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(elPart, gridSize));

                obj.solverContext_.Gxx = magPartXFFT + (1 - forwardShiftX) .* (1 - backwardShiftX) .* elPartFFT;

                obj.solverContext_.detG = obj.solverContext_.Gxx;
            elseif all(isnan(magPartX), 'all') && ~all(isnan(magPartY), 'all')
                obj.solverContext_.rooftopContent = 'onlyY';

                magPartYFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartY, gridSize));
                elPartFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(elPart, gridSize));

                obj.solverContext_.Gyy = magPartYFFT + (1 - forwardShiftY) .* (1 - backwardShiftY) .* elPartFFT;

                obj.solverContext_.detG = obj.solverContext_.Gyy;
            else
                obj.solverContext_.rooftopContent = 'all';

                magPartXFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartX, gridSize));
                magPartYFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartY, gridSize));
                elPartFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(elPart, gridSize));

                obj.solverContext_.Gxx = magPartXFFT + (1 - forwardShiftX) .* (1 - backwardShiftX) .* elPartFFT;
                obj.solverContext_.Gxy = (1 - forwardShiftX) .* (1 - backwardShiftY) .* elPartFFT;
                obj.solverContext_.Gyx = (1 - forwardShiftY) .* (1 - backwardShiftX) .* elPartFFT;
                obj.solverContext_.Gyy = magPartYFFT + (1 - forwardShiftY) .* (1 - backwardShiftY) .* elPartFFT;

                obj.solverContext_.detG = obj.solverContext_.Gxx .* obj.solverContext_.Gyy - obj.solverContext_.Gxy .* obj.solverContext_.Gyx;
            end
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentCG(obj, excitationVectorX, excitationVectorY, tolerance, maxIterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            [currentVector, converged, residueNorm, iterationCount] = obj.conjugateGradientSolver(excitationVector, tolerance, maxIterationCount);
            %[currentVector, converged, residueNorm, iterationCount] = bicgstab(@(excitationVector) obj.applyOperatorsVector(excitationVector), excitationVector, tolerance, maxIterationCount);
            %[currentVector, converged, residueNorm, iterationCount] = cgs(@(excitationVector) obj.applyOperatorsVector(excitationVector), excitationVector, tolerance, maxIterationCount);

            currentVectorX = reshape(currentVector(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentPCG(obj, excitationVectorX, excitationVectorY, tolerance, maxIterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            %[currentVector, converged, residueNorm, iterationCount] = bicgstab(@(excitationVector) obj.applyOperatorsVector(excitationVector), excitationVector, tolerance, maxIterationCount, @(r) obj.applyPreconditionerVector(r));
            [currentVector, converged, residueNorm, iterationCount] = obj.preconditionedConjugateGradientSolver(excitationVector, tolerance, maxIterationCount);

            currentVectorX = reshape(currentVector(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentBiCGSTAB(obj, excitationVectorX, excitationVectorY, tolerance, maxIterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            [currentVector, converged, residueNorm, iterationCount] = bicgstab(@(excitationVector) obj.applyOperatorsVector(excitationVector), excitationVector, tolerance, maxIterationCount);

            currentVectorX = reshape(currentVector(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentPGMRES(obj, excitationVectorX, excitationVectorY, restart, tolerance, maxIterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            [currentVector, converged, residueNorm, iterationCount] = gmres(@(excitationVector) obj.applyOperatorsVector(excitationVector), excitationVector, restart, tolerance, maxIterationCount, @(r) obj.applyPreconditionerVector(r));

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

                % just for debugging purposes
                if mod(column, 100) == 0
                    %fprintf("Computed %d columns!\n", column);
                end
            end
        end

        function preconditionerMatrix = computePreconditionerMatrix(obj)
            basisFunctionCount = obj.solverContext_.basisCountX.x * obj.solverContext_.basisCountX.y + obj.solverContext_.basisCountY.x * obj.solverContext_.basisCountY.y;
            preconditionerMatrix = zeros(basisFunctionCount);

            for column = 1:basisFunctionCount
                canonicalVector = zeros(basisFunctionCount, 1);
                canonicalVector(column) = 1;

                preconditionerMatrix(:, column) = obj.applyPreconditionerVector(canonicalVector);

                % just for debugging purposes
                if mod(column, 100) == 0
                    %fprintf("Computed %d columns!\n", column);
                end
            end
        end

        % TODO: currently not working for very big impedance matrices, use computeImpedanceMatrix1
        function impedanceMatrix = computeImpedanceMatrix2(obj)
            basisFunctionCount = obj.solverContext_.basisCountX.x * obj.solverContext_.basisCountX.y + obj.solverContext_.basisCountY.x * obj.solverContext_.basisCountY.y;

            I = eye(basisFunctionCount);

            impedanceMatrix = obj.applyOperatorsMatrix(I);
        end
    end

    methods (Access = private)
        function result = applyOperatorsVector(obj, currentVector)
            switch obj.solverContext_.rooftopContent
                case 'onlyX'
                    currentVectorX = currentVector;

                    currentVectorXFFT = fft2(currentVectorX, obj.solverContext_.paddingX, obj.solverContext_.paddingY);

                    GxFFT = obj.solverContext_.Gxx .* currentVectorXFFT;

                    Gx = ifft2(GxFFT);

                    GxSection = Gx(1:obj.solverContext_.basisCountX.x, 1:obj.solverContext_.basisCountX.y);

                    result = [GxSection(:)];
                case 'onlyY'
                    currentVectorY = currentVector.';

                    currentVectorYFFT = fft2(currentVectorY, obj.solverContext_.paddingX, obj.solverContext_.paddingY);

                    GyFFT = obj.solverContext_.Gyy .* currentVectorYFFT;

                    Gy = ifft2(GyFFT);

                    GySection = Gy(1:obj.solverContext_.basisCountY.x, 1:obj.solverContext_.basisCountY.y);

                    result = [GySection(:)];
                case 'all'
                    currentVectorX = reshape(currentVector(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
                    currentVectorY = reshape(currentVector((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);

                    paddingX = obj.solverContext_.paddingX;
                    paddingY = obj.solverContext_.paddingY;

                    currentVectorXFFT = fft2(currentVectorX, paddingX, paddingY);
                    currentVectorYFFT = fft2(currentVectorY, paddingX, paddingY);

                    GxFFT = obj.solverContext_.Gxx .* currentVectorXFFT + obj.solverContext_.Gxy .* currentVectorYFFT; 
                    GyFFT = obj.solverContext_.Gyx .* currentVectorXFFT + obj.solverContext_.Gyy .* currentVectorYFFT; 

                    Gx = ifft2(GxFFT);
                    Gy = ifft2(GyFFT);

                    GxSection = Gx(1:obj.solverContext_.basisCountX.x, 1:obj.solverContext_.basisCountX.y);
                    GySection = Gy(1:obj.solverContext_.basisCountY.x, 1:obj.solverContext_.basisCountY.y);

                    result = [GxSection(:); GySection(:)];
            end
        end

        function z = applyPreconditionerVector(obj, r)
            switch obj.solverContext_.rooftopContent
                case 'onlyX'
                    rX = r;

                    rXFFT = fft2(rX, obj.solverContext_.paddingX, obj.solverContext_.paddingY);

                    zXFFT = rXFFT ./ obj.solverContext_.detG;

                    zX = ifft2(zXFFT);

                    zXSection = zX(1:obj.solverContext_.basisCountX.x, 1:obj.solverContext_.basisCountX.y);

                    z = [zXSection(:)];
                case 'onlyY'
                    rY = r.';

                    rYFFT = fft2(rY, obj.solverContext_.paddingX, obj.solverContext_.paddingY);

                    zYFFT = rYFFT ./ obj.solverContext_.detG;

                    zY = ifft2(zYFFT);

                    zYSection = zY(1:obj.solverContext_.basisCountY.x, 1:obj.solverContext_.basisCountY.y);

                    z = [zYSection(:)];
                case 'all'
                    rX = reshape(r(1:obj.solverContext_.offsetX), obj.solverContext_.basisCountX.x, obj.solverContext_.basisCountX.y);
                    rY = reshape(r((obj.solverContext_.offsetX + 1):end), obj.solverContext_.basisCountY.x, obj.solverContext_.basisCountY.y);

                    rXFFT = fft2(rX, obj.solverContext_.paddingX, obj.solverContext_.paddingY);
                    rYFFT = fft2(rY, obj.solverContext_.paddingX, obj.solverContext_.paddingY);

                    zXFFT = (obj.solverContext_.Gyy .* rXFFT - obj.solverContext_.Gxy .* rYFFT) ./ obj.solverContext_.detG;
                    zYFFT = (obj.solverContext_.Gxx .* rYFFT - obj.solverContext_.Gyx .* rXFFT) ./ obj.solverContext_.detG;

                    zX = ifft2(zXFFT);
                    zY = ifft2(zYFFT);

                    zXSection = zX(1:obj.solverContext_.basisCountX.x, 1:obj.solverContext_.basisCountX.y);
                    zYSection = zY(1:obj.solverContext_.basisCountY.x, 1:obj.solverContext_.basisCountY.y);

                    z = [zXSection(:); zYSection(:)];
            end
        end

        function result = applyOperatorsMatrix(obj, currentMatrix)
            % Unpack parameters
            bcx = obj.solverContext_.basisCountX;
            bcy = obj.solverContext_.basisCountY;
            padX = obj.solverContext_.paddingX;
            padY = obj.solverContext_.paddingY;
            totalX = bcx.x * bcx.y;
            totalY = bcy.x * bcy.y;
            
            % TODO: implement as a general setting
            batchSize = 1024;

            switch obj.solverContext_.rooftopContent
                case 'onlyX'
                    N = size(currentMatrix, 2);
                    result = zeros(totalX, N); % Preallocate output
                    
                    % Process in batches
                    for start = 1:batchSize:N
                        %fprintf("batch processing\n")
                        endIdx = min(start + batchSize - 1, N);
                        batch = start:endIdx;
                        currentBatch = currentMatrix(:, batch);
                        
                        % Reshape and compute FFT with padding
                        X_batch = reshape(currentBatch, bcx.x, bcx.y, []);
                        X_FFT = fft2(X_batch, padX, padY); % [padX x padY x batchSize]
                        
                        % Apply Green's function and invert FFT
                        GxFFT = obj.solverContext_.Gxx .* X_FFT;
                        Gx = ifft2(GxFFT);
                        
                        % Extract and store
                        GxSection = Gx(1:bcx.x, 1:bcx.y, :);
                        result(:, batch) = reshape(GxSection, totalX, []);
                    end

                case 'onlyY'
                    N = size(currentMatrix, 2);
                    result = zeros(totalY, N); % Preallocate output
                    
                    % Process in batches
                    for start = 1:batchSize:N
                        %fprintf("batch processing\n")
                        endIdx = min(start + batchSize - 1, N);
                        batch = start:endIdx;
                        currentBatch = currentMatrix(:, batch);
                        
                        % Reshape and compute FFT with padding
                        Y_batch = reshape(currentBatch, bcy.x, bcy.y, []);
                        Y_FFT = fft2(Y_batch, padX, padY); 
                        
                        % Apply Green's function and invert FFT
                        GyFFT = obj.solverContext_.Gyy .* Y_FFT;
                        Gy = ifft2(GyFFT);
                        
                        % Extract and store
                        GySection = Gy(1:bcy.x, 1:bcy.y, :);
                        result(:, batch) = reshape(GySection, totalY, []);
                    end

                otherwise % 'all'
                    N = size(currentMatrix, 2);
                    result = zeros(totalX + totalY, N); 
                    
                    for start = 1:batchSize:N
                        endIdx = min(start + batchSize - 1, N);
                        batch = start:endIdx;
                        currentBatch = currentMatrix(:, batch);
                        
                        % Split into X/Y components
                        X_batch = reshape(currentBatch(1:totalX, :), bcx.x, bcx.y, []);
                        Y_batch = reshape(currentBatch(totalX+1:end, :), bcy.x, bcy.y, []);
                        
                        % Compute FFTs
                        X_FFT = fft2(X_batch, padX, padY);
                        Y_FFT = fft2(Y_batch, padX, padY);
                        
                        % Apply Green's functions
                        GxFFT = obj.solverContext_.Gxx .* X_FFT + obj.solverContext_.Gxy .* Y_FFT;
                        GyFFT = obj.solverContext_.Gyx .* X_FFT + obj.solverContext_.Gyy .* Y_FFT;
                        
                        % Invert FFTs and store
                        Gx = ifft2(GxFFT);
                        Gy = ifft2(GyFFT);
                        result(1:totalX, batch) = reshape(Gx(1:bcx.x, 1:bcx.y, :), totalX, []);
                        result(totalX+1:end, batch) = reshape(Gy(1:bcy.x, 1:bcy.y, :), totalY, []);
                    end
            end
        end

        function [currentVector, converged, residueNorm, iterationCount] = conjugateGradientSolver(obj, excitationVector, tolerance, maxIterationCount)
            currentVector = obj.applyPreconditionerVector(excitationVector);

            residue = excitationVector - obj.applyOperatorsVector(currentVector);
            direction = residue;
            rTr = residue.' * residue;
            excitationNorm = norm(excitationVector);
            converged = false;

            for iterationCount = 1:maxIterationCount
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

        function [currentVector, converged, residueNorm, iterationCount] = preconditionedConjugateGradientSolver(obj, excitationVector, tolerance, maxIterationCount)
            currentVector = zeros(size(excitationVector));

            residue = excitationVector - obj.applyOperatorsVector(currentVector);
            preconditionedResidue = obj.applyPreconditionerVector(residue);
            direction = preconditionedResidue;
            rTz = residue.' * preconditionedResidue;
            excitationNorm = norm(excitationVector);
            converged = false;

            for iterationCount = 1:maxIterationCount
                GDirection = obj.applyOperatorsVector(direction);
                alpha = rTz / (direction.' * GDirection);
                currentVector = currentVector + alpha * direction;
                residue = residue - alpha * GDirection;
                preconditionedResidue = obj.applyPreconditionerVector(residue);
                rTzNew = residue.' * preconditionedResidue;
                residueNorm = sqrt(rTzNew) / excitationNorm;

                if residueNorm < tolerance
                    converged = true;
                    break;
                end

                beta = rTzNew / rTz;
                direction = preconditionedResidue + beta * direction;
                rTz = rTzNew;
            end
        end
    end
end
