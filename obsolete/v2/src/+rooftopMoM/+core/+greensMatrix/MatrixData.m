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
            if nargin < 2
                options = struct();
            end

            if ~isfield(options, 'kernelSizeX') && ~isfield(options, 'kernelSizeY')
                [options.kernelSizeX, options.kernelSizeY] = obj.minimalKernelSize();
            end

            [parts, options] = obj.prepareParts_(options);

            circulants = cell(1, numel(parts));
            for n = 1:numel(parts)
                circulants{n} = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(parts{n}, options.kernelSizeX, options.kernelSizeY);
            end

            spectralMatrixData = rooftopMoM.core.greensMatrix.utils.generateSpectralGreensMatrix( ...
                circulants{1}, circulants{2}, circulants{3}, obj.rooftopContent_, options.faceSizeX, options.faceSizeY ...
            );
        end

        function spectralMatrixData = generateSplitSpectralGreensMatrix(obj, options)
        % generateSplitSpectralGreensMatrix: Generates the subband spectra for the SplitFFTAcceleratedApplicator.
        %
        %   Input:
        %       options: Optional struct with the fields modifier (applied to every part), faceSizeX and faceSizeY.
        %
        %   Output:
        %       spectralMatrixData: 1xS array of SpectralMatrixData (see utils.generateSplitSpectralGreensMatrix).

            if nargin < 2
                options = struct();
            end

            [parts, options] = obj.prepareParts_(options);

            spectralMatrixData = rooftopMoM.core.greensMatrix.utils.generateSplitSpectralGreensMatrix( ...
                parts{1}, parts{2}, parts{3}, obj.rooftopContent_, options.faceSizeX, options.faceSizeY ...
            );
        end
    end

    methods (Access = private)
        function [parts, options] = prepareParts_(obj, options)
        % prepareParts_: Applies the optional modifier to the Toeplitz generators and defaults the face sizes.
        %
        %   Input:
        %       options: Struct with the optional fields modifier, faceSizeX and faceSizeY.
        %
        %   Outputs:
        %       parts: 1x3 cell of the (modified) generators {magneticXX, magneticYY, electric}.
        %       options: The input options with both face sizes set to 1 if neither was given.

            if ~isfield(options, 'faceSizeX') && ~isfield(options, 'faceSizeY')
                options.faceSizeX = 1;
                options.faceSizeY = 1;
            end

            parts = {obj.magneticXX_, obj.magneticYY_, obj.electric_};

            if isfield(options, 'modifier')
                parts = cellfun(options.modifier, parts, 'UniformOutput', false);
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
