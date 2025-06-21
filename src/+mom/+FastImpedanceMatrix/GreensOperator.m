classdef GreensOperator < handle
    properties (Access = protected)
        strategyContext_;

        greensOperator_;
    end

    methods
        function obj = ConjugateGradientSolver(strategyContext)
            obj.strategyContext_ = strategyContext;

            %fftw('planner','exhaustive');
        end
    end

    methods (Access = public)
        function generateOperators(obj, chunkMemoryMB)
            % generate the irreducible parts of magnetic and electric part
            magPartStrategy = mom.FastImpedanceMatrix.ComputeMagPart();
            elPartStrategy = mom.FastImpedanceMatrix.ComputeElPart();

            [magPartX, magPartY] = magPartStrategy.compute(obj.strategyContext_, chunkMemoryMB);
            elPart = elPartStrategy.compute(obj.strategyContext_, chunkMemoryMB);

            gridSize = double(obj.strategyContext_.rooftopBasisFunctions.getGridSize());

            % determine the rooftop counts in each direction
            obj.greensOperator_.basisCountX.x = gridSize(1) - 1;
            obj.greensOperator_.basisCountX.y = gridSize(2);
            obj.greensOperator_.basisCountY.x = gridSize(1);
            obj.greensOperator_.basisCountY.y = gridSize(2) - 1;
            obj.greensOperator_.offsetX = (gridSize(1) - 1) * gridSize(2);

            M = gridSize(1) * 2;
            N = gridSize(2) * 2;

            obj.greensOperator_.paddingX = M;
            obj.greensOperator_.paddingY = N;

            % compute shift operators
            [u, v] = ndgrid(0:(M - 1), 0:(N - 1));
            forwardShiftX = exp((1i * pi * u) / gridSize(1));
            forwardShiftY = exp((1i * pi * v) / gridSize(2));
            backwardShiftX = conj(forwardShiftX);
            backwardShiftY = conj(forwardShiftY);

            % construct the greens operator is spectra
            if ~all(isnan(magPartX), 'all') && all(isnan(magPartY), 'all')
                obj.greensOperator_.rooftopContent = 'onlyX';

                magPartXFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartX, gridSize));
                elPartFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(elPart, gridSize));

                obj.greensOperator_.Gxx = magPartXFFT + (1 - forwardShiftX) .* (1 - backwardShiftX) .* elPartFFT;

                obj.greensOperator_.detG = obj.greensOperator_.Gxx;
            elseif all(isnan(magPartX), 'all') && ~all(isnan(magPartY), 'all')
                obj.greensOperator_.rooftopContent = 'onlyY';

                magPartYFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartY, gridSize));
                elPartFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(elPart, gridSize));

                obj.greensOperator_.Gyy = magPartYFFT + (1 - forwardShiftY) .* (1 - backwardShiftY) .* elPartFFT;

                obj.greensOperator_.detG = obj.greensOperator_.Gyy;
            else
                obj.greensOperator_.rooftopContent = 'all';

                magPartXFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartX, gridSize));
                magPartYFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(magPartY, gridSize));
                elPartFFT = fft2(mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(elPart, gridSize));

                obj.greensOperator_.Gxx = magPartXFFT + (1 - forwardShiftX) .* (1 - backwardShiftX) .* elPartFFT;
                obj.greensOperator_.Gxy = (1 - forwardShiftX) .* (1 - backwardShiftY) .* elPartFFT;
                obj.greensOperator_.Gyx = (1 - forwardShiftY) .* (1 - backwardShiftX) .* elPartFFT;
                obj.greensOperator_.Gyy = magPartYFFT + (1 - forwardShiftY) .* (1 - backwardShiftY) .* elPartFFT;

                obj.greensOperator_.detG = obj.greensOperator_.Gxx .* obj.greensOperator_.Gyy - obj.greensOperator_.Gxy .* obj.greensOperator_.Gyx;
            end

            obj.greensOperator_.invDetG = 1 ./ obj.greensOperator_.detG;
        end
    end

    methods (Access = public)
        function solverContext = getSolverContext(obj)
            solverContext = obj.greensOperator_;
        end
    end
end
