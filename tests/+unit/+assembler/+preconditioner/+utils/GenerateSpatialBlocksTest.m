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
            generators = rooftopMoM.assembler.greensMatrix.FullAssembler().generate(context);
            matrix = testCase.generateDenseMatrix(generators, rooftops);

            blocks = cell(2);
            [blocks{1, 1}, blocks{1, 2}, blocks{2, 1}, blocks{2, 2}] = ...
                rooftopMoM.assembler.preconditioner.utils.generateSpatialBlocks(generators, 0.1, 0.15);

            directions = 'xy';
            counts = [rooftops.countDir('x'), rooftops.countDir('y')];
            first = [0, counts(1)];

            for a = find(counts)
                for b = find(counts)
                    expected = matrix(first(a) + (1:counts(a)), first(b) + (1:counts(b)));
                    actual = assemble_(blocks{a, b}, rooftops.faceOffsets(directions(a), 1:counts(a)), ...
                        rooftops.faceOffsets(directions(b), 1:counts(b)), grid);

                    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-13 * max(abs(matrix(:))));
                end
            end
        end
    end
end

function block = assemble_(generator, rowOffsets, columnOffsets, grid)
    % The entry of two rooftops is the generator at the offset between their first faces.
    p = rowOffsets(1, :).' - columnOffsets(1, :);
    q = rowOffsets(2, :).' - columnOffsets(2, :);
    block = generator(sub2ind(size(generator), p + grid(1), q + grid(2)));
end
