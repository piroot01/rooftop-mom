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
            % assembleImpl_: Assemble both magnetic xx/yy parts and electric part,
            %   embed them to circulant tensors, take FFT and combine with
            %   the appropriate shift operators to produce Fourier-domain kernels.

            primitive = context.rooftops().primitive();
            [gridSizeX, gridSizeY] = primitive.gridSize();

            magneticAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.chunkPlanner);
            electricAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.chunkPlanner);

            % Compute the magnetic/electric part toeplitz kernels.
            [magneticXX, magneticYY] = magneticAssembler.generateFull(context);
            electric = electricAssembler.generateFull(context);

            % Embed Toeplitz kernels into circulant tensors and compute FFTs.
            gmxxCircFFT = obj.embedToeplitzAndFourier_(magneticXX, gridSizeX, gridSizeY);
            gmyyCircFFT = obj.embedToeplitzAndFourier_(magneticYY, gridSizeX, gridSizeY);
            geCircFFT = obj.embedToeplitzAndFourier_(electric, gridSizeX, gridSizeY);

            % Compute forward/backward shift operators for both axes.
            [sizeX, sizeY] = size(geCircFFT);
            [forwardShiftX, backwardShiftX, forwardShiftY, backwardShiftY] = ...
                obj.computeShiftOperators_(sizeX, sizeY);

            % Combine the parts into Fourier-domain kernels.
            [gxxFFT, gxyFFT, gyxFFT, gyyFFT] = ...
                obj.combineParts_(gmxxCircFFT, gmyyCircFFT, geCircFFT, ...
                    forwardShiftX, backwardShiftX, ...
                    forwardShiftY, backwardShiftY);

            output = rooftopMoM.core.greensMatrix.MatrixData( ...
                'gxxFFT', gxxFFT, ...
                'gxyFFT', gxyFFT, ...
                'gyxFFT', gyxFFT, ...
                'gyyFFT', gyyFFT ...
            );
        end

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
            sizeX = size(geCircFFT, 1);
            [forwardShiftX, backwardShiftX, ~, ~] = obj.computeShiftOperators_(sizeX, 1);

            % Build Gxx Fourier kernel.
            gxxFFT = gmxxCircFFT + (1 - forwardShiftX) .* (1 - backwardShiftX) .* geCircFFT;

            output = rooftopMoM.core.greensMatrix.MatrixData( ...
                'gxxFFT', gxxFFT ...
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
            geCircFFT   = obj.embedToeplitzAndFourier_(electric,   gridSizeX, gridSizeY);

            % Only Y-direction shifts are needed.
            sizeY = size(geCircFFT, 2);
            [~, ~, forwardShiftY, backwardShiftY] = obj.computeShiftOperators_(1, sizeY);

            % Build Gyy Fourier kernel.
            gyyFFT = gmyyCircFFT + (1 - forwardShiftY) .* (1 - backwardShiftY) .* geCircFFT;

            output = rooftopMoM.core.greensMatrix.MatrixData( ...
                'gyyFFT', gyyFFT ...
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

        function [forwardX, backwardX, forwardY, backwardY] = computeShiftOperators_(~, sizeX, sizeY)
            % computeShiftOperators_: Compute phase multipliers corresponding to spatial shifts in
            %   the Fourier domain. The shift used here is exp(i*pi*u/sizeX)
            %   and similarly for Y. The backward shift is the complex
            %   conjugate of the forward shift.

            % If one of the sizes equals 1 we still create ndgrid for
            % consistent broadcasting with elementwise ops.
            [u, v] = ndgrid(0:(sizeX - 1), 0:(sizeY - 1));

            forwardX = exp((1i * pi * u) / sizeX);
            forwardY = exp((1i * pi * v) / sizeY);

            backwardX = conj(forwardX);
            backwardY = conj(forwardY);
        end

        function [gxx, gxy, gyx, gyy] = combineParts_(~, gmxxFFT, gmyyFFT, geFFT, fShiftX, bShiftX, fShiftY, bShiftY)
            % combinePartsAll_: Combine magnetic and electric Fourier-domain contributions
            %   into the full set of Fourier-domain Green operator: Gxx, Gxy, Gyx, Gyy.
            %
            %    Formula (elementwise):
            %       Gxx = Gm_xx + (1 - fX).*(1 - bX).*Ge
            %       Gxy =            (1 - fX).*(1 - bY).*Ge
            %       Gyx =            (1 - fY).*(1 - bX).*Ge
            %       Gyy = Gm_yy + (1 - fY).*(1 - bY).*Ge
            %
            %   NOTE: All inputs are already FFT(embedded-circulant) tensors.

            gxx = gmxxFFT + (1 - fShiftX) .* (1 - bShiftX) .* geFFT;
            gxy = (1 - fShiftX) .* (1 - bShiftY) .* geFFT;
            gyx = (1 - fShiftY) .* (1 - bShiftX) .* geFFT;
            gyy = gmyyFFT + (1 - fShiftY) .* (1 - bShiftY) .* geFFT;
        end
    end
end
