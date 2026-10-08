classdef MeshUniformRooftopsGenome < rooftopMoM.basisFunctions.MeshUniformRooftops
    properties (Access = private)
        genome_ logical;
    end

    methods
        function obj = MeshUniformRooftopsGenome(mesh, quadrature, genome)
            obj@rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quadrature);

            obj.genome_ = genome;
        end
    end

    methods
        function generate(obj)
            generate@rooftopMoM.basisFunctions.MeshDivergenceConforming(obj);

            obj.filterIndices_();

            obj.quadrature().generate();
        end
    end

    methods (Access = private)
        function filterIndices_(obj)
            % NOTE: We assume that the genome respect the nature that first are the rooftops in the X direction and then in the Y direction.
            % NOTE: This rooftops works only on a single mesh not a mesh cluster yet, this might be a good thing to add in the future.
            indicesX = obj.indicesX_{1};
            indicesY = obj.indicesY_{1};

            nX = obj.countDir('x');
            nY = obj.countDir('y');

            genomeX = obj.genome_(1:nX);
            genomeY = obj.genome_((nX + 1):(nX + nY));

            obj.indicesX_{1} = indicesX(:, genomeX);
            obj.indicesY_{1} = indicesY(:, genomeY);
        end
    end
end
