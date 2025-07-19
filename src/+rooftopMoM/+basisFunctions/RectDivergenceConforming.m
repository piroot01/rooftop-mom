classdef RectDivergenceConforming < rooftopMoM.basisFunctions.BasisFunctionsInterface
    properties (Access = private)
        mesh_;

        indicesX_;
        indicesY_;
    end

    methods
        function obj = RectDivergenceConforming(mesh)
            obj.mesh_ = mesh;
        end
    end

    methods
        function generate(obj, genome)
            indices = rooftopMoM.mesh.utils.generateMeshInteriorEdges(obj.mesh_);

            % Apply the genome mask if available.
            if ~(nargin < 2 || isempty(genome))
                indices = cellfun(@(indices, genome) indices(:, genome(:)), indices, genome, 'UniformOutput', false);
            end

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
    end
end
