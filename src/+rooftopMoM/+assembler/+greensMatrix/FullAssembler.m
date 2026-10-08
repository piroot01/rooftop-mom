classdef FullAssembler < rooftopMoM.assembler.greensMatrix.Assembler
% FullAssembler: Assembles the Green's matrix generators with the tensor Gauss-Legendre rule on every face pair.
%   The static part 1 / R is integrated analytically on coincident faces.

    methods
        function generators = generate(~, context)
        % generate: Returns the generators of the rooftop directions present on the grid.
        %
        %   Input:
        %       context: rooftopMoM.assembler.AssemblerContext.
        %
        %   Output:
        %       generators: rooftopMoM.assembler.greensMatrix.Generators.

            start = tic;
            rooftops = context.rooftops();
            magneticXX = [];
            magneticYY = [];

            if rooftops.countDir('x') > 0
                magneticXX = rooftopMoM.assembler.greensMatrix.generateMagneticPart(context, 'x');
            end

            if rooftops.countDir('y') > 0
                magneticYY = rooftopMoM.assembler.greensMatrix.generateMagneticPart(context, 'y');
            end

            electric = rooftopMoM.assembler.greensMatrix.generateElectricPart(context);

            generators = rooftopMoM.assembler.greensMatrix.Generators(magneticXX, magneticYY, electric, ...
                rooftops.content());

            rooftopMoM.logging.Logger.instance().info('Assembled the generators of %d rooftops in %.3f s (%.1f MB).', ...
                rooftops.count(), toc(start), generators.memory() / 2^20);
        end
    end
end
