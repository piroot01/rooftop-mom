classdef EstimateMemoryTest < unit.Base
% EstimateMemoryTest: The memory estimates of rooftopMoM.resources against the objects they describe, once built.

    properties (TestParameter)
        % countX, countY: plates, strips, the smallest plate and sizes that need a larger fast FFT.
        grid = {[4, 3], [7, 1], [1, 6], [2, 2], [9, 5], [13, 11]};
        product = {'fft', 'splitFFT', 'leanSplitFFT'};
    end

    methods (Test)
        function testGenerators(testCase, grid)
            context = testCase.createContext(grid(1), grid(2), 0.1, 0.15, 2, 2 * pi);
            generators = rooftopMoM.assembler.greensMatrix.MomentAssembler().generate(context);

            testCase.verifyEqual(rooftopMoM.resources.estimateGenerators(context.rooftops()), generators.memory());
        end

        function testProduct(testCase, grid, product)
            [rooftops, generators] = testCase.assemble_(grid);
            applicator = testCase.createProduct(applicatorName_(product), rooftops, generators);

            [stored, working] = rooftopMoM.resources.estimateProduct(rooftops, product);

            testCase.verifyEqual([stored, working], [applicator.memory(), applicator.workingMemory()]);
        end

        function testPreconditioner(testCase, grid)
            [rooftops, generators] = testCase.assemble_(grid);
            preconditioner = rooftopMoM.assembler.preconditioner.generateCirculant(rooftops, generators, ...
                @rooftopMoM.assembler.preconditioner.kernels.generateTChan, ...
                @rooftopMoM.assembler.preconditioner.generatePartwiseSpectra);

            [stored, working] = rooftopMoM.resources.estimatePreconditioner(rooftops);

            testCase.verifyEqual([stored, working], [preconditioner.memory(), preconditioner.workingMemory()]);
        end
    end

    methods (Access = private)
        function [rooftops, generators] = assemble_(testCase, grid)
            context = testCase.createContext(grid(1), grid(2), 0.1, 0.15, 2, 2 * pi);
            rooftops = context.rooftops();
            generators = rooftopMoM.assembler.greensMatrix.MomentAssembler().generate(context);
        end
    end
end

function name = applicatorName_(product)
    names = struct('fft', 'FFTApplicator', 'splitFFT', 'SplitFFTApplicator', 'leanSplitFFT', ...
        'LeanSplitFFTApplicator');
    name = names.(product);
end
