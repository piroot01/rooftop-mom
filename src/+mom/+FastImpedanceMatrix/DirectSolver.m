classdef DirectSolver < mom.FastImpedanceMatrix.GreensOperator
    methods (Access = public)
        function obj = DirectSolver(strategyContext)
            obj = obj@mom.FastImpedanceMatrix.GreensOperator(strategyContext);
        end
    end

    methods (Access = public)
        function currentVector = solve(obj, excitationVector)
            excitationVectorX = reshape(excitationVector(1:obj.greensOperator_.offsetX), obj.greensOperator_.basisCountX.x, obj.greensOperator_.basisCountX.y);
            excitationVectorY = reshape(excitationVector((obj.greensOperator_.offsetX + 1):end), obj.greensOperator_.basisCountY.x, obj.greensOperator_.basisCountY.y);

            excitationVectorXFFT = fft2(excitationVectorX, obj.greensOperator_.paddingX, obj.greensOperator_.paddingY);
            excitationVectorYFFT = fft2(excitationVectorY, obj.greensOperator_.paddingX, obj.greensOperator_.paddingY);

            currentVectorXFFT = (obj.greensOperator_.Gyy .* excitationVectorXFFT - obj.greensOperator_.Gxy .* excitationVectorYFFT) .* obj.greensOperator_.invDetG;
            currentVectorYFFT = (obj.greensOperator_.Gxx .* excitationVectorYFFT - obj.greensOperator_.Gyx .* excitationVectorXFFT) .* obj.greensOperator_.invDetG;

            currentVectorX = ifft2(currentVectorXFFT);
            currentVectorY = ifft2(currentVectorYFFT);

            currentVectorXSection = currentVectorX(1:obj.greensOperator_.basisCountX.x, 1:obj.greensOperator_.basisCountX.y);
            currentVectorYSection = currentVectorY(1:obj.greensOperator_.basisCountY.x, 1:obj.greensOperator_.basisCountY.y);

            currentVector = [currentVectorXSection(:); currentVectorYSection(:)];
        end
    end
end
