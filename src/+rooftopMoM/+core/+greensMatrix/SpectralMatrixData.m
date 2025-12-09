classdef SpectralMatrixData < handle
    properties (Access = private)
        blockXX_;
        blockXY_;
        blockYX_;
        blockYY_;
    end

    methods
        function obj = SpectralMatrixData(varargin)
            p = inputParser();

            addParameter(p, 'blockXX', []);
            addParameter(p, 'blockXY', []);
            addParameter(p, 'blockYX', []);
            addParameter(p, 'blockYY', []);

            parse(p, varargin{:});

            obj.blockXX_ = p.Results.blockXX;
            obj.blockXY_ = p.Results.blockXY;
            obj.blockYX_ = p.Results.blockYX;
            obj.blockYY_ = p.Results.blockYY;
        end
    end

    methods
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
