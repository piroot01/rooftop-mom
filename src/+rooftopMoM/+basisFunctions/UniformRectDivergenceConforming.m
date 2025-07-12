classdef UniformRectDivergenceConforming < handle
    properties (Access = private)
        mesh_;

        indicesX_;
        indicesY_;
    end

    methods
        function obj = UniformRectDivergenceConforming(mesh)
            obj.mesh_ = mesh;
        end
    end

    methods
        function generate(obj)
            indices = rooftopMoM.mesh.utils.generateMeshInteriorEdges(obj.mesh_);

            obj.indicesX_ = cellfun(@(xs) rooftopMoM.mesh.utils.filterInteriorEdgesByDir(xs, 'x'), indices, 'UniformOutput', false);
            obj.indicesY_ = cellfun(@(xs) rooftopMoM.mesh.utils.filterInteriorEdgesByDir(xs, 'y'), indices, 'UniformOutput', false);
        end

        function count = count(obj)
            count = sum(cellfun(@(xs) size(xs, 2), obj.indicesX_)) + ...
                    sum(cellfun(@(xs) size(xs, 2), obj.indicesY_));
        end

        function count = countDir(obj, dir)
            if dir == 'x'
                count = sum(cellfun(@(xs) size(xs, 2), obj.indicesX_));
            elseif dir == 'y'
                count = sum(cellfun(@(xs) size(xs, 2), obj.indicesY_));
            else
                count = obj.count();
            end
        end

        function indices = indices(obj)
            indicesCell = obj.indicesCell();
            indices = [indicesCell{:}];
        end

        function indices = indicesCell(obj, shifted)
            if nargin < 2 || isempty(shifted)
                shifted = true;
            end

            meshCount = obj.mesh_.count();

            indicesX = obj.indicesX_;
            indicesY = obj.indicesY_;

            % If shifted==false, this is a zero vector; otherwise mesh offsets.
            offsets = uint32(shifted) * obj.mesh_.offset();

            indices = arrayfun(@(n) [indicesX{n}, indicesY{n}] + offsets(n), 1:meshCount, 'UniformOutput', false);
        end

        function indices = indicesCellDir(obj, dir, shifted)
            if nargin < 3 || isempty(shifted)
                shifted = true;
            end

            meshCount = obj.mesh_.count();

            % If shifted==false, this is a zero vector; otherwise mesh offsets.
            offsets = uint32(shifted) * obj.mesh_.offset();

            if dir == 'x'
                indices = obj.indicesX_;
            elseif dir == 'y'
                indices = obj.indicesY_;
            else
                indices = obj.indicesCell(false);
            end

            indices = arrayfun(@(n) indices{n} + offsets(n), 1:meshCount, 'UniformOutput', false);
        end

        function mesh = mesh(obj)
            mesh = obj.mesh_;
        end
    end
end
