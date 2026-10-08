classdef Spectra
% Spectra: The 2x2 block spectra of the Green's matrix on a frequency grid. XX and YY act within the x- and the
%   y-directed rooftops, XY and YX couple them. The blocks of a missing rooftop direction are empty.

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
