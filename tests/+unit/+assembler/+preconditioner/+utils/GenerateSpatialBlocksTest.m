classdef GenerateSpatialBlocksTest < unit.Base
% GenerateSpatialBlocksTest: Tests of rooftopMoM.assembler.preconditioner.utils.generateSpatialBlocks.

    properties (TestParameter)
        % countX, countY: a plate, a grid two faces wide and the strips along x and y.
        grid = {[4, 3], [2, 3], [5, 1], [1, 4]};
    end

    methods (Test)
        function testBlocksReproduceTheGreensMatrix(testCase, grid)
            context = testCase.createContext(grid(1), grid(2), 0.1, 0.15, 3, 2 * pi);
            rooftops = context.rooftops();
            generators = rooftopMoM.assembler.greensMatrix.MomentAssembler().generate(context);
            matrix = testCase.generateDenseMatrix(generators, rooftops);

            blocks = cell(2);
            [blocks{1, 1}, blocks{1, 2}, blocks{2, 1}, blocks{2, 2}] = ...
                rooftopMoM.assembler.preconditioner.utils.generateSpatialBlocks(generators, 0.1, 0.15);

            directions = 'xy';
            counts = [rooftops.countDir('x'), rooftops.countDir('y')];
            first = [0, counts(1)];

            for m = find(counts)
                for n = find(counts)
                    expected = matrix(first(m) + (1:counts(m)), first(n) + (1:counts(n)));
                    actual = assemble_(blocks{m, n}, rooftops.faceOffsets(directions(m), 1:counts(m)), ...
                        rooftops.faceOffsets(directions(n), 1:counts(n)), grid);

                    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-13 * max(abs(matrix(:))));
                end
            end
        end
    end
end

function block = assemble_(generator, rowOffsets, columnOffsets, grid)
    % The entry of two rooftops is the generator at the offset between their first faces.
    m = rowOffsets(1, :).' - columnOffsets(1, :);
    n = rowOffsets(2, :).' - columnOffsets(2, :);
    block = generator(sub2ind(size(generator), m + grid(1), n + grid(2)));
end
