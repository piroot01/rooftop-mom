classdef FullAssembler < rooftopMoM.core.greensMatrix.Assembler
% FullAssembler: Assembles the full Green's operator in the Fourier (circulant) domain
%   for rooftop basis functions. This class decides which directional
%   parts (xx, yy or both) are required based on the primitive grid and
%   then orchestrates assembly by delegating to magnetic and electric
%   part assemblers and a small set of private helper methods.

    properties (Access = private)
        chunkPlanner  % planner object used by the part assemblers
    end

    methods
        function obj = FullAssembler(chunkPlanner)
            % constructor: store planner (keeps class small)
            obj.chunkPlanner = chunkPlanner;
        end

        function output = generateFull(obj, context)
            % generateFull: Decide whether to assemble both directions
            %   or only one, according to primitive grid size.

            primitive = context.rooftops().primitive();
            [gridSizeX, gridSizeY] = primitive.gridSize();

            if gridSizeX ~= 1 && gridSizeY ~= 1
                % Both directions present -> full assembly.
                output = obj.generateImpl_(context);
                return;
            elseif gridSizeX ~= 1 && gridSizeY == 1
                % Only X direction nontrivial -> assemble XX part only.
                output = obj.generateXXImpl_(context);
                return;
            elseif gridSizeX == 1 && gridSizeY ~= 1
                % Only Y direction nontrivial -> assemble YY part only.
                output = obj.generateYYImpl_(context);
                return;
            else
                % Both are 1 (degenerate) -> error (illegal primitive).
                error('The dimensions of the primitive are illegal %dx%d.', gridSizeX, gridSizeY);
            end
        end

        function generateSpecified(~, ~)
            error('generateSpecified is not implemented yet.');
        end
    end
    
    methods (Access = private)
        function output = generateImpl_(obj, context)
            % generateImpl_: Assemble both magnetic xx/yy parts and electric part,

            magneticAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.chunkPlanner);
            electricAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.chunkPlanner);

            % Compute the magnetic/electric part toeplitz kernels.
            [magneticXX, magneticYY] = magneticAssembler.generateFull(context);
            electric = electricAssembler.generateFull(context);

            output = rooftopMoM.core.greensMatrix.MatrixData( ...
                'magneticXX', magneticXX, ...
                'magneticYY', magneticYY, ...
                'electric', electric, ...
                'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.BOTH ...
            );
        end
        %{
            % the output here will be the MatrixData instance, the instance should contain the rooftopCOntent info
            % in the GreensMatrixApplicator in the initialization (generate method) there will be:
            %   (i) FFTAcceleratedApplicator:
            %       - function that takes the rooftopContent and matrixData returns the potentially three embedded circulant matrices
            %         if the matrix is not required then it will be empty or NaN ???
            %       - question is whether the next function should take the embedded or plus fft2 (for now we will include in the subsequent func)
            %       - then there will be the function that inputs the potentially three embedded matrices and constructs the potentially
            %         four parts of the spectrall greensMatrix

            % Embed Toeplitz kernels into circulant tensors and compute FFTs.
            gmxxCircFFT = obj.embedToeplitzAndFourier_(magneticXX, gridSizeX, gridSizeY);
            gmyyCircFFT = obj.embedToeplitzAndFourier_(magneticYY, gridSizeX, gridSizeY);
            geCircFFT = obj.embedToeplitzAndFourier_(electric, gridSizeX, gridSizeY);

            % Compute forward/backward shift operators for both axes.
            [forwardShiftX, backwardShiftX, forwardShiftY, backwardShiftY] = ...
                obj.computeShiftOperators_(gridSizeX - 1, gridSizeY - 1);

            % Combine the parts into Fourier-domain kernels.
            [faceSizeX, faceSizeY] = primitive.faceSize();
            gxxFFT = gmxxCircFFT + (1 / faceSizeX)^2 .* (1 - forwardShiftX) .* (backwardShiftX - 1) .* geCircFFT;
            gxyFFT = (1 / (faceSizeX * faceSizeY)) .* (1 - forwardShiftX) .* (backwardShiftY - 1) .* geCircFFT;
            gyxFFT = (1 / (faceSizeY * faceSizeX)) .* (1 - forwardShiftY) .* (backwardShiftX - 1) .* geCircFFT;
            gyyFFT = gmyyCircFFT + (1 / faceSizeY)^2 .* (1 - forwardShiftY) .* (backwardShiftY - 1) .* geCircFFT;

            determinant = gxxFFT .* gyyFFT - gxyFFT .* gyxFFT;

            output = rooftopMoM.core.greensMatrix.MatrixData( ...
                'gxxFFT', gxxFFT, ...
                'gxyFFT', gxyFFT, ...
                'gyxFFT', gyxFFT, ...
                'gyyFFT', gyyFFT, ...
                'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.BOTH, ...
                'paddingX', 2 * (gridSizeX - 1), ...
                'paddingY', 2 * (gridSizeY - 1), ...
                'determinant', determinant ...
            );
        end
        %}

        function output = generateXXImpl_(obj, context)
            % generateXXImpl_: Assemble only the XX directional operator (magnetic XX + electric).

            primitive = context.rooftops().primitive();
            [gridSizeX, gridSizeY] = primitive.gridSize();

            magneticAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.chunkPlanner);
            electricAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.chunkPlanner);

            magneticXX = magneticAssembler.generateFullDir(context, 'x');
            electric = electricAssembler.generateFull(context);

            % embed + FFT
            gmxxCircFFT = obj.embedToeplitzAndFourier_(magneticXX, gridSizeX, gridSizeY);
            geCircFFT = obj.embedToeplitzAndFourier_(electric, gridSizeX, gridSizeY);

            % Only X-direction shifts are needed.
            [forwardShiftX, backwardShiftX, ~, ~] = obj.computeShiftOperators_(gridSizeX - 1, 1);

            % Build Gxx Fourier kernel.
            faceSizeX = primitive.sizeX();
            gxxFFT = gmxxCircFFT + (1 / faceSizeX)^2 .* (1 - forwardShiftX) .* (backwardShiftX - 1) .* geCircFFT;

            % The determinant is just gxxFFT so we wont store it once more.

            output = rooftopMoM.core.greensMatrix.MatrixData( ...
                'gxxFFT', gxxFFT, ...
                'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.XONLY, ... 
                'paddingX', 2 * (gridSizeX - 1), ...
                'paddingY', 2 ...
            );
        end

        function output = generateYYImpl_(obj, context)
            % generateYYImpl_: Assemble only the YY directional operator (magnetic YY + electric).

            primitive = context.rooftops().primitive();
            [gridSizeX, gridSizeY] = primitive.gridSize();

            magneticAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.chunkPlanner);
            electricAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.chunkPlanner);

            magneticYY = magneticAssembler.generateFullDir(context, 'y');
            electric = electricAssembler.generateFull(context);

            % embed + FFT
            gmyyCircFFT = obj.embedToeplitzAndFourier_(magneticYY, gridSizeX, gridSizeY);
            geCircFFT = obj.embedToeplitzAndFourier_(electric, gridSizeX, gridSizeY);

            % Only Y-direction shifts are needed.
            [~, ~, forwardShiftY, backwardShiftY] = obj.computeShiftOperators_(1, gridSizeY - 1);

            % Build Gyy Fourier kernel.
            faceSizeY = primitive.sizeY();
            gyyFFT = gmyyCircFFT + (1 / faceSizeY)^2 .* (1 - forwardShiftY) .* (backwardShiftY - 1) .* geCircFFT;

            % The determinant is just gxxFFT so we wont store it once more.

            output = rooftopMoM.core.greensMatrix.MatrixData( ...
                'gyyFFT', gyyFFT, ...
                'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.YONLY, ... 
                'paddingX', 2, ...
                'paddingY', 2 * (gridSizeY - 1) ...
            );
        end
    end

    methods (Access = private)
        function circFFT = embedToeplitzAndFourier_(~, toeplitzKernel, gridSizeX, gridSizeY)
            % embedToeplitzAndFourier_: Embeds a Toeplitz kernel into a circulant tensor and returns its 2D Fourier transform (fft2).

            circulantTensor = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant( ...
                toeplitzKernel, gridSizeX, gridSizeY ...
            );

            circFFT = fft2(circulantTensor);
        end

        function [forwardX, backwardX, forwardY, backwardY] = computeShiftOperators_(~, M, N)
            % computeShiftOperators_: Compute phase multipliers corresponding to spatial shifts in
            %   the Fourier domain. The shift used here is exp(i*pi*u/sizeX)
            %   and similarly for Y. The backward shift is the complex
            %   conjugate of the forward shift.

            M = double(M);
            N = double(N);

            [u, v] = ndgrid(0:(2 * M - 1), 0:(2 * N - 1));

            forwardX = exp((1i * pi * u) / M);
            forwardY = exp((1i * pi * v) / N);

            backwardX = conj(forwardX);
            backwardY = conj(forwardY);
        end
    end
end
