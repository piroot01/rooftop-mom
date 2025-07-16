classdef MagneticPartAssembler < rooftopMoM.core.impedanceMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_ rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner;
    end

    methods
        function obj = MagneticPartAssembler(assemblerChunkPlanner)
            obj.assemblerChunkPlanner_ = assemblerChunkPlanner;
        end
    end

    methods
        function Z = generate(obj, context)
            chunks = obj.assemblerChunkPlanner_.generateFull(context, true);

            Z = obj.generateChunks(context, chunks);
        end

        function Z = generateChunks(obj, context, chunks)
            wavenumber = context.wavenumber();
            matrixCount = numel(wavenumber);
            basisSize = context.rooftops().count();

            % Generate the quadrature, rooftop, rooftopDiv points (same for every matrix).
            [points, values] = context.rooftops().generateRooftops();
            [quadraturePoints, rooftopPoints] = rooftopMoM.core.impedanceMatrix.utils.convertRooftopPointsToArrays(points);

            sin

            % Pre-allocate the final impedance matrix (matrices).
            Z = complex(zeros(basisSize, basisSize, matrixCount));

            for n = 1:matrixCount
                Z(:, :, n) = obj.generateImpl_(context, chunks, wavenumber(n));
            end
        end
    end

    methods (Access = private)
        function Z = generateImpl_(obj, context, wavenumber)
            p1 = [1, 1, 2, 2];
            p2 = [1, 2, 1, 2];
            pCount = numel(p1);
        end
    end
end
