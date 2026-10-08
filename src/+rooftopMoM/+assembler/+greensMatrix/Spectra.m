classdef Spectra
% Spectra: The spectra of a 2x2 block circulant matrix on a frequency grid, such as the circulant embedding of the
%   Green's matrix. XX and YY act within the x- and the y-directed rooftops, XY and YX couple them. The blocks of a
%   missing rooftop direction are empty.

    properties (Access = private)
        blockXX_;
        blockXY_;
        blockYX_;
        blockYY_;
    end

    methods
        function obj = Spectra(blockXX, blockXY, blockYX, blockYY)
            obj.blockXX_ = blockXX;
            obj.blockXY_ = blockXY;
            obj.blockYX_ = blockYX;
            obj.blockYY_ = blockYY;
        end
    end

    methods
        function determinant = determinant(obj, content)
        % determinant: Returns the determinant of the 2x2 blocks at every frequency, or the single block of a strip.
        %
        %   Input:
        %       content: rooftopMoM.basisFunctions.RooftopContent.
        %
        %   Output:
        %       determinant: Array of the size of the frequency grid.

            if content == rooftopMoM.basisFunctions.RooftopContent.XONLY
                determinant = obj.blockXX_;
            elseif content == rooftopMoM.basisFunctions.RooftopContent.YONLY
                determinant = obj.blockYY_;
            else
                determinant = obj.blockXX_ .* obj.blockYY_ - obj.blockXY_ .* obj.blockYX_;
            end
        end

        function inverse = invert(obj, content)
        % invert: Returns the spectra of the inverse of the block circulant matrix. The FFT diagonalizes every
        %   circulant block, so the inverse is the inverse of the 2x2 matrix [XX, XY; YX, YY] at every frequency, the
        %   reciprocal of the single block for a strip.
        %
        %   Input:
        %       content: rooftopMoM.basisFunctions.RooftopContent.
        %
        %   Output:
        %       inverse: rooftopMoM.assembler.greensMatrix.Spectra.

            if content == rooftopMoM.basisFunctions.RooftopContent.XONLY
                inverse = rooftopMoM.assembler.greensMatrix.Spectra(1 ./ obj.blockXX_, [], [], []);
            elseif content == rooftopMoM.basisFunctions.RooftopContent.YONLY
                inverse = rooftopMoM.assembler.greensMatrix.Spectra([], [], [], 1 ./ obj.blockYY_);
            else
                determinant = obj.determinant(content);
                inverse = rooftopMoM.assembler.greensMatrix.Spectra(obj.blockYY_ ./ determinant, ...
                    -obj.blockXY_ ./ determinant, -obj.blockYX_ ./ determinant, obj.blockXX_ ./ determinant);
            end
        end

        function [sizeX, sizeY] = gridSize(obj)
        % gridSize: Returns the size of the frequency grid.

            blocks = {obj.blockXX_, obj.blockYY_, obj.blockXY_, obj.blockYX_};
            [sizeX, sizeY] = size(blocks{find(~cellfun(@isempty, blocks), 1)});
        end

        function bytes = memory(obj)
        % memory: Returns the bytes of the blocks of all elements of the (array) object.

            bytes = 0;

            for n = 1:numel(obj)
                bytes = bytes + rooftopMoM.assembler.utils.generateByteCount(obj(n).blockXX_, obj(n).blockXY_, ...
                    obj(n).blockYX_, obj(n).blockYY_);
            end
        end

        function blockXX = blockXX(obj)
            blockXX = obj.blockXX_;
        end

        function blockXY = blockXY(obj)
            blockXY = obj.blockXY_;
        end

        function blockYX = blockYX(obj)
            blockYX = obj.blockYX_;
        end

        function blockYY = blockYY(obj)
            blockYY = obj.blockYY_;
        end
    end
end
