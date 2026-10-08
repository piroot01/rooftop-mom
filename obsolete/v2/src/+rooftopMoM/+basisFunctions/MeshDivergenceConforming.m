classdef MeshDivergenceConforming < rooftopMoM.basisFunctions.MeshBasisFunctionsInterface
    properties (Access = protected)
        mesh_;

        indicesX_;
        indicesY_;
    end

    methods
        function obj = MeshDivergenceConforming(mesh)
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

        function counts = counts(obj)
            counts = cellfun(@(xs, ys) size(xs, 2) + size(ys, 2), obj.indicesX_, obj.indicesY_);
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
            meshCount = obj.mesh_.count(); 

            indicesX = obj.indicesX_;
            indicesY = obj.indicesY_;

            indices = arrayfun(@(n) [indicesX{n}, indicesY{n}], 1:meshCount, 'UniformOutput', false);
        end

        function indices = indicesDir(obj, dir)
            meshCount = obj.mesh_.count();

            if dir == 'x'
                indices = obj.indicesX_;
            elseif dir == 'y'
                indices = obj.indicesY_;
            else
                indices = obj.indices();
            end

            indices = arrayfun(@(n) indices{n}, 1:meshCount, 'UniformOutput', false);
        end

        function mesh = mesh(obj)
            mesh = obj.mesh_;
        end

        function memory = memory(obj)
            indicesXCount = sum(cellfun(@numel, obj.indicesX_));
            indicesYCount = sum(cellfun(@numel, obj.indicesY_));

            memory = (indicesXCount + indicesYCount) * rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.UInt32);
        end
    end
end
