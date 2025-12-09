classdef MatrixData < handle
    properties (Access = private)
        magneticXX_;
        magneticYY_;
        electric_;

        rooftopContent_ rooftopMoM.basisFunctions.RooftopContent;
    end

    methods
        function obj = MatrixData(varargin)
            p = inputParser();

            addParameter(p, 'magneticXX', []);
            addParameter(p, 'magneticYY', []);
            addParameter(p, 'electric', []);

            addParameter(p, 'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.NONE);

            parse(p, varargin{:});

            obj.magneticXX_ = p.Results.magneticXX;
            obj.magneticYY_ = p.Results.magneticYY;
            obj.electric_ = p.Results.electric;

            obj.rooftopContent_ = p.Results.rooftopContent;
        end
    end

    methods
        function spectralMatrixData = generateSpectralGreensMatrix(obj, options)
            if nargin < 2 || (~isfield(options, 'kernelSizeX') && ~isfield(options, 'kernelSizeY'))
                % Size Extraction
                if ~isempty(obj.electric_)
                    [kernelSizeX, kernelSizeY] = size(obj.electric_);
                elseif ~isempty(obj.magneticXX_)
                    [kernelSizeX, kernelSizeY] = size(obj.magneticXX_);
                elseif ~isempty(obj.magneticYY_)
                    [kernelSizeX, kernelSizeY] = size(obj.magneticYY_);
                else
                    error('At least one of the kernels must be nonempty!');
                end

                options.kernelSizeX = kernelSizeX;
                options.kernelSizeY = kernelSizeY;
            end

            if nargin < 2 || (~isfield(options, 'faceSizeX') && ~isfield(options, 'faceSizeY'))
                options.faceSizeX = 1;
                options.faceSizeY = 1;
            end

            [magneticXXCirc, magneticYYCirc, electricCirc] = obj.generateCirculants_(options.kernelSizeX, options.kernelSizeY);

            [blockXX, blockXY, blockYX, blockYY] = obj.generateSpectralGreensMatrixBlocks_( ...
                magneticXXCirc, magneticYYCirc, electricCirc, options.faceSizeX, options.faceSizeY ...
            );

            spectralMatrixData = rooftopMoM.core.greensMatrix.SpectralMatrixData( ...
                'blockXX', blockXX, ...
                'blockXY', blockXY, ...
                'blockYX', blockYX, ...
                'blockYY', blockYY ...
            );
        end
    end

    methods (Access = private)
        function [magneticXXCirc, magneticYYCirc, electricCirc] = generateCirculants_(obj, kernelSizeX, kernelSizeY)
            magneticXXCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(obj.magneticXX_, kernelSizeX, kernelSizeY);
            magneticYYCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(obj.magneticYY_, kernelSizeX, kernelSizeY);
            electricCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(obj.electric_, kernelSizeX, kernelSizeY);
        end

        function [blockXX, blockXY, blockYX, blockYY] = generateSpectralGreensMatrixBlocks_(obj, magneticXXCirc, magneticYYCirc, electricCirc, faceSizeX, faceSizeY)
            % Size Extraction
            if ~isempty(electricCirc)
                [circSizeX, circSizeY] = size(electricCirc);
            elseif ~isempty(magneticXXCirc)
                [circSizeX, circSizeY] = size(magneticXXCirc);
            elseif ~isempty(magneticYYCirc)
                [circSizeX, circSizeY] = size(magneticYYCirc);
            else
                error('At least one of the circulants must be nonempty!');
            end

            % Compute Spectra (FFT) & Sanitize
            magneticXXCircSpectr = fft2(magneticXXCirc);
            if isempty(magneticXXCircSpectr)
                magneticXXCircSpectr = 0;
            end

            magneticYYCircSpectr = fft2(magneticYYCirc);
            if isempty(magneticYYCircSpectr)
                magneticYYCircSpectr = 0;
            end

            electricCircSpectr = fft2(electricCirc);
            if isempty(electricCircSpectr)
                electricCircSpectr = 0;
            end

            % Generate Shift Operators (Calculate once, use everywhere)
            [forwardShiftX, backwardShiftX, forwardShiftY, backwardShiftY] = ...
                rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(circSizeX, circSizeY);

            % Determine Active Modes
            isXActive = (obj.rooftopContent_ == rooftopMoM.basisFunctions.RooftopContent.XONLY) || ...
                (obj.rooftopContent_ == rooftopMoM.basisFunctions.RooftopContent.BOTH);

            isYActive = (obj.rooftopContent_ == rooftopMoM.basisFunctions.RooftopContent.YONLY) || ...
                (obj.rooftopContent_ == rooftopMoM.basisFunctions.RooftopContent.BOTH);

            % Compute Blocks
            blockXX = []; blockXY = []; blockYX = []; blockYY = [];

            % Compute XX Component
            if isXActive
                termMagnetic = magneticXXCircSpectr;
                termElectric = (1 / faceSizeX)^2 .* (1 - forwardShiftX) .* (backwardShiftX - 1) .* electricCircSpectr;
                blockXX = termMagnetic + termElectric;
            end

            % Compute YY Component
            if isYActive
                termMagnetic = magneticYYCircSpectr;
                termElectric = (1 / faceSizeY)^2 .* (1 - forwardShiftY) .* (backwardShiftY - 1) .* electricCircSpectr;
                blockYY = termMagnetic + termElectric;
            end

            % Compute Cross-Coupling (XY / YX)
            if isXActive && isYActive
                blockXY = (1 / (faceSizeX * faceSizeY)) .* (1 - forwardShiftX) .* (backwardShiftY - 1) .* electricCircSpectr;
                blockYX = (1 / (faceSizeY * faceSizeX)) .* (1 - forwardShiftY) .* (backwardShiftX - 1) .* electricCircSpectr;
            end
        end
    end

    methods
        function magneticXX = magneticXX(obj)
            magneticXX = obj.magneticXX_;
        end

        function magneticYY = magneticYY(obj)
            magneticYY = obj.magneticYY_;
        end

        function electric = electric(obj)
            electric = obj.electric_;
        end

        function rooftopContent = rooftopContent(obj)
            rooftopContent = obj.rooftopContent_;
        end
    end
end
