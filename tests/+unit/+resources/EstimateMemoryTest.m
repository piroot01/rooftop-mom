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

        function testChunk(testCase, grid)
            % The estimate is the smallest chunk with which each pass of the assembly runs in one chunk.
            context = testCase.createContext(grid(1), grid(2), 0.1, 0.15, [3, 2], 2 * pi);
            bytes = rooftopMoM.resources.estimateChunk(context.rooftops(), context.quadratureOrders());

            testCase.verifyEqual(testCase.countChunks_(grid, bytes), [1, 1]);
            testCase.verifyGreaterThan(max(testCase.countChunks_(grid, bytes - 1)), 1);
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
        function counts = countChunks_(testCase, grid, chunkSize)
        % countChunks_: Returns the numbers of chunks of the near and of the far pass of an assembly, read from its
        %   trace log.

            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.TRACE, {sink});
            context = testCase.createContext(grid(1), grid(2), 0.1, 0.15, [3, 2], 2 * pi, chunkSize);
            rooftopMoM.assembler.greensMatrix.MomentAssembler().generate(context);

            totals = regexp(strjoin(sink.lines(), newline), 'Integrated chunk 1 of (\d+)\.', 'tokens');
            counts = cellfun(@(token) str2double(token{1}), totals);
        end

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
