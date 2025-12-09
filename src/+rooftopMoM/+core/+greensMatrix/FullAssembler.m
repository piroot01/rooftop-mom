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

        function matrixData = generate(obj, context)
            primitive = context.rooftops().primitive();
            [gridSizeX, gridSizeY] = primitive.gridSize();

            if gridSizeX ~= 1 && gridSizeY ~= 1
                % Both directions present -> full assembly.
                matrixData = obj.generateImpl_(context);
                return;
            elseif gridSizeX ~= 1 && gridSizeY == 1
                % Only X direction nontrivial -> assemble XX part only.
                matrixData = obj.generateXXImpl_(context);
                return;
            elseif gridSizeX == 1 && gridSizeY ~= 1
                % Only Y direction nontrivial -> assemble YY part only.
                matrixData = obj.generateYYImpl_(context);
                return;
            else
                % Both are 1 (degenerate) -> error (illegal primitive).
                error('The dimensions of the primitive are illegal %dx%d.', gridSizeX, gridSizeY);
            end
        end

        function [magneticXX, magneticYY, electric] = generateFull(obj, context)
            magneticAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.chunkPlanner);
            electricAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.chunkPlanner);

            [magneticXX, magneticYY] = magneticAssembler.generateFull(context);
            electric = electricAssembler.generateFull(context);
        end

        function generateSpecified(~, ~)
            error('generateSpecified is not implemented yet.');
        end
    end
    
    methods (Access = private)
        function matrixData = generateImpl_(obj, context)
        % generateImpl_: Assemble both magnetic xx/yy parts and electric part,

            magneticAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.chunkPlanner);
            electricAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.chunkPlanner);

            % Compute the magnetic/electric part toeplitz kernels.
            [magneticXX, magneticYY] = magneticAssembler.generateFull(context);
            electric = electricAssembler.generateFull(context);

            matrixData = rooftopMoM.core.greensMatrix.MatrixData( ...
                'magneticXX', magneticXX, ...
                'magneticYY', magneticYY, ...
                'electric', electric, ...
                'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.BOTH ...
            );
        end
        
        function matrixData = generateXXImpl_(obj, context)
        % generateXXImpl_: Assemble only the XX directional operator (magnetic XX + electric).

            magneticAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.chunkPlanner);
            electricAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.chunkPlanner);

            magneticXX = magneticAssembler.generateFullDir(context, 'x');
            electric = electricAssembler.generateFull(context);

            matrixData = rooftopMoM.core.greensMatrix.MatrixData( ...
                'magneticXX', magneticXX, ...
                'magneticYY', [], ...
                'electric', electric, ...
                'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.XONLY ...
            );
        end
        
        function matrixData = generateYYImpl_(obj, context)
            % generateYYImpl_: Assemble only the YY directional operator (magnetic YY + electric).

            magneticAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.chunkPlanner);
            electricAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.chunkPlanner);

            magneticYY = magneticAssembler.generateFullDir(context, 'y');
            electric = electricAssembler.generateFull(context);

            matrixData = rooftopMoM.core.greensMatrix.MatrixData( ...
                'magneticXX', [], ...
                'magneticYY', magneticYY, ...
                'electric', electric, ...
                'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.YONLY ...
            );
        end
    end
end
