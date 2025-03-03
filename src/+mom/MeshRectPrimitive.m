classdef MeshRectPrimitive < handle
    properties (Access = private)
        patchCountX_ uint32;
        patchCountY_ uint32;
        patchSizeX_ double;
        patchSizeY_ double;
        genom_ logical;
    end

    methods
        function obj = mesh(patchCountX, patchCountY, patchSizeX, patchSizeY, genom)
            obj.patchCountX_ = patchCountX;
            obj.patchCountY_ = patchCountY;
            obj.patchSizeX_ = patchSizeX;
            obj.patchSizeY_ = patchSizeY;

            if ~isempty(genom) && numel(genom) == getPatchCount(obj)
                obj.genom_ = logical(genom);
            else
                obj.genom_ = true(1, patch_count_x_direction * patch_count_y_direction);
            end
        end
    end

    methods
        function patchCountX = getPatchCountX(obj)
            patchCountX = obj.patchCountX_;
        end

        function patchCountY = getPatchCountY(obj)
            patchCountY = obj.patchCountY_;
        end

        function patchSizeX = getPatchSizeX(obj)
            patchSizeX = obj.patchSizeX_;
        end

        function patchSizeY =getPatchSizeY(obj)
            patchSizeY = obj.patchSizeY_;
        end

        function count = getPatchCount(obj)
            count = obj.patchCountX_ * obj.patchCountY_;
        end

        function indices = getAllIndices(obj)
            indices = 1:getPatchCount(obj);
        end

        function indices = getIndices(obj)
            allIndices = obj.getAllIndices();
            indices = allIndices(obj.genom_);
        end
    end
end