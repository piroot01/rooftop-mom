classdef PerturbedConjugateGradientSolver < mom.FastImpedanceMatrix.PerturbedGreensOperator
    methods (Access = public)
        function obj = PerturbedConjugateGradientSolver(strategyContext)
            obj = obj@mom.FastImpedanceMatrix.PerturbedGreensOperator(strategyContext);
        end
    end

    methods (Access = public)
        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentCG(obj, excitationVectorX, excitationVectorY, tolerance, maxIterationCount, initialCurrentVectorXY)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            if nargin == 5 || isempty(initialCurrentVectorXY)
                initialCurrentVectorXY_ = zeros(size(excitationVector)); 
            elseif nargin == 6
                initialCurrentVectorXY_ = initialCurrentVectorXY;
            end

            [currentVector, converged, residueNorm, iterationCount] = obj.conjugateGradientSolver(excitationVector, tolerance, maxIterationCount, initialCurrentVectorXY_);

            currentVectorX = reshape(currentVector(1:obj.greensOperator_.offsetX), obj.greensOperator_.basisCountX.x, obj.greensOperator_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.greensOperator_.offsetX + 1):end), obj.greensOperator_.basisCountY.x, obj.greensOperator_.basisCountY.y);
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentPCG(obj, excitationVectorX, excitationVectorY, tolerance, maxIterationCount, initialCurrentVectorXY)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            if nargin == 5 || isempty(initialCurrentVectorXY)
                initialCurrentVectorXY_ = zeros(size(excitationVector)); 
            elseif nargin == 6
                initialCurrentVectorXY_ = initialCurrentVectorXY;
            end

            [currentVector, converged, residueNorm, iterationCount] = obj.preconditionedConjugateGradientSolver(excitationVector, tolerance, maxIterationCount, initialCurrentVectorXY_);

            currentVectorX = reshape(currentVector(1:obj.greensOperator_.offsetX), obj.greensOperator_.basisCountX.x, obj.greensOperator_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.greensOperator_.offsetX + 1):end), obj.greensOperator_.basisCountY.x, obj.greensOperator_.basisCountY.y);
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentBiCGSTAB(obj, excitationVectorX, excitationVectorY, tolerance, maxIterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            [currentVector, converged, residueNorm, iterationCount] = bicgstab(@(excitationVector) obj.applyOperatorsVector(excitationVector), excitationVector, tolerance, maxIterationCount);

            currentVectorX = reshape(currentVector(1:obj.greensOperator_.offsetX), obj.greensOperator_.basisCountX.x, obj.greensOperator_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.greensOperator_.offsetX + 1):end), obj.greensOperator_.basisCountY.x, obj.greensOperator_.basisCountY.y);
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentPGMRES(obj, excitationVectorX, excitationVectorY, restart, tolerance, maxIterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            [currentVector, converged, residueNorm, iterationCount] = gmres(@(excitationVector) obj.applyOperatorsVector(excitationVector), excitationVector, restart, tolerance, maxIterationCount, @(r) obj.applyPreconditionerVector(r));

            currentVectorX = reshape(currentVector(1:obj.greensOperator_.offsetX), obj.greensOperator_.basisCountX.x, obj.greensOperator_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.greensOperator_.offsetX + 1):end), obj.greensOperator_.basisCountY.x, obj.greensOperator_.basisCountY.y);
        end

        function [currentVectorX, currentVectorY, converged, residueNorm, iterationCount] = solveCurrentGMRES(obj, excitationVectorX, excitationVectorY, restart, tolerance, maxIterationCount)
            excitationVector = [excitationVectorX(:); excitationVectorY(:)];

            [currentVector, converged, residueNorm, iterationCount] = gmres(@(excitationVector) obj.applyOperatorsVector(excitationVector), excitationVector, restart, tolerance, maxIterationCount);

            currentVectorX = reshape(currentVector(1:obj.greensOperator_.offsetX), obj.greensOperator_.basisCountX.x, obj.greensOperator_.basisCountX.y);
            currentVectorY = reshape(currentVector((obj.greensOperator_.offsetX + 1):end), obj.greensOperator_.basisCountY.x, obj.greensOperator_.basisCountY.y);
        end

        function impedanceMatrix = computeImpedanceMatrix1(obj)
            basisFunctionCount = obj.greensOperator_.basisCountX.x * obj.greensOperator_.basisCountX.y + obj.greensOperator_.basisCountY.x * obj.greensOperator_.basisCountY.y;
            impedanceMatrix = zeros(basisFunctionCount);

            for column = 1:basisFunctionCount
                canonicalVector = zeros(basisFunctionCount, 1);
                canonicalVector(column) = 1;

                impedanceMatrix(:, column) = obj.applyOperatorsVector(canonicalVector);

                % just for debugging purposes
                if mod(column, 100) == 0
                    fprintf("[impedance matrix]: Computed %d columns!\n", column);
                end
            end
        end

        function radiationResistanceMatrix = computeRadiationResistanceMatrix(obj)
            % we make use of the fact that the matrix is symmetric
            basisFunctionCount = obj.greensOperator_.basisCountX.x * obj.greensOperator_.basisCountX.y + obj.greensOperator_.basisCountY.x * obj.greensOperator_.basisCountY.y;
            radiationResistanceMatrix = zeros(basisFunctionCount);

            for column = 1:basisFunctionCount
                canonicalVector = zeros(basisFunctionCount, 1);
                canonicalVector(column) = 1;

                radiationResistanceMatrix(:, column) = real(obj.applyOperatorsVector(canonicalVector));

                % just for debugging purposes
                if mod(column, 100) == 0
                    fprintf("[radiationResistanceMatrix]: Computed %d columns!\n", column);
                end
            end
        end

        function radiationResistanceMatrix = computeParRadiationResistanceMatrix(obj)
            bx = obj.greensOperator_.basisCountX;
            by = obj.greensOperator_.basisCountY;
            N  = bx.x*bx.y + by.x*by.y;

            radiationResistanceMatrix = zeros(N);

            pObj = parallel.pool.Constant(obj);

            parfor column = 1:N
                localObj = pObj.Value;

                canonicalVector = zeros(N, 1);
                canonicalVector(column) = 1;

                radiationResistanceMatrix(:, column) = real(localObj.applyOperatorsVector(canonicalVector));

                % just for debugging purposes
                if mod(column,100) == 0
                    fprintf('[radiationResistanceMatrix]: Computed %d/%d columns on worker %d\n', ...
                            column, N, getCurrentTask().ID);
                end
            end
        end


        function preconditionerMatrix = computePreconditionerMatrix(obj)
            basisFunctionCount = obj.greensOperator_.basisCountX.x * obj.greensOperator_.basisCountX.y + obj.greensOperator_.basisCountY.x * obj.greensOperator_.basisCountY.y;
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

        function Pr = computeRadiatedPower(obj, I)
            v = obj.applyOperatorsVector(I);

            Pr = 1/2 * real(I' * v);
        end

        % TODO: currently not working for very big impedance matrices, use computeImpedanceMatrix1
        function impedanceMatrix = computeImpedanceMatrix2(obj)
            basisFunctionCount = obj.greensOperator_.basisCountX.x * obj.greensOperator_.basisCountX.y + obj.greensOperator_.basisCountY.x * obj.greensOperator_.basisCountY.y;

            I = eye(basisFunctionCount);

            impedanceMatrix = obj.applyOperatorsMatrix(I);
        end

        function [totalBytes, totalMiB] = getOperatorSizeMiB(obj)
            Gfields = {'Gxx','Gxy','Gyx','Gyy'};
            
            totalBytes = 0;
            for k = 1:numel(Gfields)
                f = Gfields{k};
                if isfield(obj.greensOperator_, f)
                    tmp = obj.greensOperator_.(f);    %#ok<NASGU>
                    S = whos('tmp');
                    totalBytes = totalBytes + S.bytes;
                end
            end
            
            % Convert to MiB (1 MiB = 2^20 bytes)
            totalMiB = totalBytes / 2^20;
        end
    end

    methods (Access = private)
        function result = applyOperatorsVector(obj, currentVector)
            switch obj.greensOperator_.rooftopContent
                case 'onlyX'
                    currentVectorX = currentVector;

                    currentVectorXFFT = fft2(currentVectorX, obj.greensOperator_.paddingX, obj.greensOperator_.paddingY);

                    GxFFT = obj.greensOperator_.Gxx .* currentVectorXFFT;

                    Gx = ifft2(GxFFT);

                    GxSection = Gx(1:obj.greensOperator_.basisCountX.x, 1:obj.greensOperator_.basisCountX.y);

                    result = [GxSection(:)];
                case 'onlyY'
                    currentVectorY = currentVector.';

                    currentVectorYFFT = fft2(currentVectorY, obj.greensOperator_.paddingX, obj.greensOperator_.paddingY);

                    GyFFT = obj.greensOperator_.Gyy .* currentVectorYFFT;

                    Gy = ifft2(GyFFT);

                    GySection = Gy(1:obj.greensOperator_.basisCountY.x, 1:obj.greensOperator_.basisCountY.y);

                    result = [GySection(:)];
                case 'all'
                    currentVectorX = reshape(currentVector(1:obj.greensOperator_.offsetX), obj.greensOperator_.basisCountX.x, obj.greensOperator_.basisCountX.y);
                    currentVectorY = reshape(currentVector((obj.greensOperator_.offsetX + 1):end), obj.greensOperator_.basisCountY.x, obj.greensOperator_.basisCountY.y);

                    paddingX = obj.greensOperator_.paddingX;
                    paddingY = obj.greensOperator_.paddingY;

                    currentVectorXFFT = fft2(currentVectorX, paddingX, paddingY);
                    currentVectorYFFT = fft2(currentVectorY, paddingX, paddingY);

                    GxFFT = obj.greensOperator_.Gxx .* currentVectorXFFT + obj.greensOperator_.Gxy .* currentVectorYFFT; 
                    GyFFT = obj.greensOperator_.Gyx .* currentVectorXFFT + obj.greensOperator_.Gyy .* currentVectorYFFT; 

                    Gx = ifft2(GxFFT);
                    Gy = ifft2(GyFFT);

                    GxSection = Gx(1:obj.greensOperator_.basisCountX.x, 1:obj.greensOperator_.basisCountX.y);
                    GySection = Gy(1:obj.greensOperator_.basisCountY.x, 1:obj.greensOperator_.basisCountY.y);

                    result = [GxSection(:); GySection(:)];
            end
        end

        function z = applyPreconditionerVector(obj, r)
            switch obj.greensOperator_.rooftopContent
                case 'onlyX'
                    rX = r;

                    rXFFT = fft2(rX, obj.greensOperator_.paddingX, obj.greensOperator_.paddingY);

                    zXFFT = rXFFT .* obj.greensOperator_.invDetG;

                    zX = ifft2(zXFFT);

                    zXSection = zX(1:obj.greensOperator_.basisCountX.x, 1:obj.greensOperator_.basisCountX.y);

                    z = [zXSection(:)];
                case 'onlyY'
                    rY = r.';

                    rYFFT = fft2(rY, obj.greensOperator_.paddingX, obj.greensOperator_.paddingY);

                    zYFFT = rYFFT .* obj.greensOperator_.invDetG;

                    zY = ifft2(zYFFT);

                    zYSection = zY(1:obj.greensOperator_.basisCountY.x, 1:obj.greensOperator_.basisCountY.y);

                    z = [zYSection(:)];
                case 'all'
                    rX = reshape(r(1:obj.greensOperator_.offsetX), obj.greensOperator_.basisCountX.x, obj.greensOperator_.basisCountX.y);
                    rY = reshape(r((obj.greensOperator_.offsetX + 1):end), obj.greensOperator_.basisCountY.x, obj.greensOperator_.basisCountY.y);

                    rXFFT = fft2(rX, obj.greensOperator_.paddingX, obj.greensOperator_.paddingY);
                    rYFFT = fft2(rY, obj.greensOperator_.paddingX, obj.greensOperator_.paddingY);

                    zXFFT = (obj.greensOperator_.Gyy .* rXFFT - obj.greensOperator_.Gxy .* rYFFT) .* obj.greensOperator_.invDetG;
                    zYFFT = (obj.greensOperator_.Gxx .* rYFFT - obj.greensOperator_.Gyx .* rXFFT) .* obj.greensOperator_.invDetG;

                    zX = ifft2(zXFFT);
                    zY = ifft2(zYFFT);

                    zXSection = zX(1:obj.greensOperator_.basisCountX.x, 1:obj.greensOperator_.basisCountX.y);
                    zYSection = zY(1:obj.greensOperator_.basisCountY.x, 1:obj.greensOperator_.basisCountY.y);

                    z = [zXSection(:); zYSection(:)];
            end
        end

        function result = applyOperatorsMatrix(obj, currentMatrix)
            % Unpack parameters
            bcx = obj.greensOperator_.basisCountX;
            bcy = obj.greensOperator_.basisCountY;
            padX = obj.greensOperator_.paddingX;
            padY = obj.greensOperator_.paddingY;
            totalX = bcx.x * bcx.y;
            totalY = bcy.x * bcy.y;
            
            % TODO: implement as a general setting
            batchSize = 1024;

            switch obj.greensOperator_.rooftopContent
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
                        GxFFT = obj.greensOperator_.Gxx .* X_FFT;
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
                        GyFFT = obj.greensOperator_.Gyy .* Y_FFT;
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
                        GxFFT = obj.greensOperator_.Gxx .* X_FFT + obj.greensOperator_.Gxy .* Y_FFT;
                        GyFFT = obj.greensOperator_.Gyx .* X_FFT + obj.greensOperator_.Gyy .* Y_FFT;
                        
                        % Invert FFTs and store
                        Gx = ifft2(GxFFT);
                        Gy = ifft2(GyFFT);
                        result(1:totalX, batch) = reshape(Gx(1:bcx.x, 1:bcx.y, :), totalX, []);
                        result(totalX+1:end, batch) = reshape(Gy(1:bcy.x, 1:bcy.y, :), totalY, []);
                    end
            end
        end

        function [currentVector, converged, residueNorm, iterationCount] = conjugateGradientSolver(obj, excitationVector, tolerance, maxIterationCount, initialCurrentVector)
            currentVector = initialCurrentVector;
            residue = excitationVector - obj.applyOperatorsVector(initialCurrentVector);
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

                if mod(iterationCount, 1) == 0
                    fprintf("Elapsed %d iterations, residueNorm = %f.\n", iterationCount, residueNorm);
                end
            end
        end

        function [currentVector, converged, residueNorm, iterationCount] = preconditionedConjugateGradientSolver(obj, excitationVector, tolerance, maxIterationCount, initialCurrentVector)
            currentVector = initialCurrentVector;
            residue = excitationVector - obj.applyOperatorsVector(initialCurrentVector);
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

                if mod(iterationCount, 1) == 0
                    fprintf("Elapsed %d iterations, residueNorm = %f.\n", iterationCount, residueNorm);
                end
            end
        end
    end
end

