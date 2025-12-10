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
        function [kernelSizeX, kernelSizeY] = minimalKernelSize(obj)
            if ~isempty(obj.electric_)
                [kernelSizeX, kernelSizeY] = size(obj.electric_);
            else
                switch obj.rooftopContent_
                    case rooftopMoM.basisFunctions.RooftopContent.XONLY
                        kernelSizeX = size(obj.magneticXX_, 1);
                        kernelSizeY = 1;
                    case rooftopMoM.basisFunctions.RooftopContent.YONLY
                        kernelSizeX = 1;
                        kernelSizeY = size(obj.magneticYY_, 2);
                    case rooftopMoM.basisFunctions.RooftopContent.BOTH
                        kernelSizeX = size(obj.magneticXX_, 2);
                        kernelSizeY = size(obj.magneticYY_, 1);
                    otherwise
                        error('Unknown rooftop content: %s', rooftopContent);
                end
            end
        end

        function spectralMatrixData = generateSpectralGreensMatrix(obj, options)
            if nargin < 2 || (~isfield(options, 'kernelSizeX') && ~isfield(options, 'kernelSizeY'))
                [options.kernelSizeX, options.kernelSizeY] = obj.minimalKernelSize();
            end

            if nargin < 2 || (~isfield(options, 'faceSizeX') && ~isfield(options, 'faceSizeY'))
                options.faceSizeX = 1;
                options.faceSizeY = 1;
            end

            [magneticXXCirc, magneticYYCirc, electricCirc] = obj.generateCirculants_(options.kernelSizeX, options.kernelSizeY);

            spectralMatrixData = rooftopMoM.core.greensMatrix.utils.generateSpectralGreensMatrix( ...
                magneticXXCirc, magneticYYCirc, electricCirc, obj.rooftopContent_, options.faceSizeX, options.faceSizeY ...
            );
        end
    end

    methods (Access = private)
        function [magneticXXCirc, magneticYYCirc, electricCirc] = generateCirculants_(obj, kernelSizeX, kernelSizeY)
            magneticXXCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(obj.magneticXX_, kernelSizeX, kernelSizeY);
            magneticYYCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(obj.magneticYY_, kernelSizeX, kernelSizeY);
            electricCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(obj.electric_, kernelSizeX, kernelSizeY);
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
