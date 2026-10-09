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
        determinant = determinant(obj, content);

        inverse = invert(obj, content);

        [sizeX, sizeY] = gridSize(obj);

        bytes = memory(obj);

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
