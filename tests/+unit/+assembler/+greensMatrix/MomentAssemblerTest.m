classdef MomentAssemblerTest < unit.Base
% MomentAssemblerTest: Tests of rooftopMoM.assembler.greensMatrix.MomentAssembler.

    methods (Test)
        function testGeneratorSizes(testCase)
            generators = testCase.generate_(5, 4, 0.1, 0.1, 3);

            testCase.verifySize(generators.magneticXX(), [4, 4]);
            testCase.verifySize(generators.magneticYY(), [5, 3]);
            testCase.verifySize(generators.electric(), [5, 4]);
            testCase.verifyEqual(generators.content(), rooftopMoM.basisFunctions.RooftopContent.BOTH);
        end

        function testStripHasOneDirection(testCase)
            generators = testCase.generate_(6, 1, 0.1, 0.1, 3);

            testCase.verifySize(generators.magneticXX(), [5, 1]);
            testCase.verifyEmpty(generators.magneticYY());
            testCase.verifyEqual(generators.content(), rooftopMoM.basisFunctions.RooftopContent.XONLY);
        end

        function testEntriesDoNotDependOnTheGridSize(testCase)
            narrow = testCase.generate_(2, 3, 0.1, 0.1, 3);
            wide = testCase.generate_(5, 3, 0.1, 0.1, 3);
            wideXX = wide.magneticXX();
            wideYY = wide.magneticYY();
            wideElectric = wide.electric();

            testCase.verifyEqual(narrow.magneticXX(), wideXX(1, :), 'RelTol', 1e-13);
            testCase.verifyEqual(narrow.magneticYY(), wideYY(1:2, :), 'RelTol', 1e-13);
            testCase.verifyEqual(narrow.electric(), wideElectric(1:2, :), 'RelTol', 1e-13);
        end

        function testTransposedGridTransposesTheGenerators(testCase)
            % Exchanging the directions maps the x-directed rooftops onto the y-directed ones.
            generators = testCase.generate_(5, 4, 0.1, 0.15, 4);
            transposed = testCase.generate_(4, 5, 0.15, 0.1, 4);

            testCase.verifyEqual(transposed.magneticYY(), generators.magneticXX().', 'RelTol', 1e-13);
            testCase.verifyEqual(transposed.electric(), generators.electric().', 'RelTol', 1e-13);
        end

        function testChunkingDoesNotChangeTheResult(testCase)
            oneChunk = testCase.createContext(5, 4, 0.1, 0.15, [6, 3], 2 * pi, 2^30);
            oneOffsetPerChunk = testCase.createContext(5, 4, 0.1, 0.15, [6, 3], 2 * pi, 1);

            whole = rooftopMoM.assembler.greensMatrix.MomentAssembler().generate(oneChunk);
            chunked = rooftopMoM.assembler.greensMatrix.MomentAssembler().generate(oneOffsetPerChunk);

            testCase.verifyEqual(chunked.magneticXX(), whole.magneticXX(), 'RelTol', 1e-14);
            testCase.verifyEqual(chunked.magneticYY(), whole.magneticYY(), 'RelTol', 1e-14);
            testCase.verifyEqual(chunked.electric(), whole.electric(), 'RelTol', 1e-14);
        end

        function testLatticeEqualsTheIntegrationOfEveryOffset(testCase)
            % Square cells and one order for both passes make the near and the far rule the same, so integrating every
            % offset on its own with that rule must give the generators that the far pass forms on the lattice.
            [counts, cellSize, order, k] = deal([9, 7], 0.1, 6, 2 * pi);
            generators = testCase.generate_(counts(1), counts(2), cellSize, cellSize, order);

            utils = 'rooftopMoM.assembler.greensMatrix.utils.';
            rule = feval([utils, 'generateReducedRule'], order, [1, 1]);
            [dX, dY] = ndgrid(0:(counts(1) - 1), 0:(counts(2) - 1));
            offsets = [dX(:).'; dY(:).'];
            scales = 1i * rooftopMoM.assembler.Constants.Z0 / (4 * pi) * cellSize^4 * [1 / k, k];
            parts = struct('electric', zeros(counts, 'like', 1i), 'magneticXX', zeros(counts - [1, 0], 'like', 1i), ...
                'magneticYY', zeros(counts - [0, 1], 'like', 1i));
            parts = feval([utils, 'accumulateMoments'], parts, ...
                feval([utils, 'integrateMoments'], rule, offsets, [cellSize, cellSize], k), offsets, scales);
            parts = feval([utils, 'accumulateMoments'], parts, ...
                rooftopMoM.assembler.singularIntegrals.generateExtractedMoments(cellSize, cellSize, k), ...
                [0, 1, 0, 1; 0, 0, 1, 1], scales);

            testCase.verifyEqual(generators.electric(), parts.electric, 'RelTol', 1e-14);
            testCase.verifyEqual(generators.magneticXX(), parts.magneticXX, 'RelTol', 1e-14);
            testCase.verifyEqual(generators.magneticYY(), parts.magneticYY, 'RelTol', 1e-14);
        end

        function testDefaultOrdersReachTheRoundingLevel(testCase)
            % Square cells, elongated cells and the cells of a thin strip: the default orders against much higher ones.
            cases = {12, 12, 0.1, 0.1; 6, 40, 0.1, 0.01; 30, 1, 0.025, 0.004};

            for n = 1:size(cases, 1)
                default = testCase.generate_(cases{n, :}, [16, 8]);
                precise = testCase.generate_(cases{n, :}, [32, 24]);

                testCase.verifyEqual(default.magneticXX(), precise.magneticXX(), 'RelTol', 3e-14);
                testCase.verifyEqual(default.magneticYY(), precise.magneticYY(), 'RelTol', 3e-14);
                testCase.verifyEqual(default.electric(), precise.electric(), 'RelTol', 3e-14);
            end
        end

        function testTouchingEntriesMatchAdaptiveIntegration(testCase)
            % The moments of the Green's function itself, integrated adaptively with integral2 in the two-dimensional
            % representation of the near-field note, without any extraction.
            [sizes, k] = deal([0.1, 0.15], 2 * pi);
            generators = testCase.generate_(3, 3, sizes(1), sizes(2), [16, 8]);
            scale = 1i * rooftopMoM.assembler.Constants.Z0 / (4 * pi) * prod(sizes)^2;
            moment = @(offset, kinds) integrateAdaptively_(sizes, k, offset, kinds);
            electric = generators.electric();
            magneticXX = generators.magneticXX();

            for offset = [0, 1, 0, 1; 0, 0, 1, 1]
                testCase.verifyEqual(electric(offset(1) + 1, offset(2) + 1), scale / k * moment(offset, [0, 0]), ...
                    'RelTol', 1e-10);
            end

            % The rooftop at the origin: the four wedge pairs meet in the cell pairs at the offsets (0, 0) and (1, 0).
            coincident = cellfun(@(kinds) moment([0, 0], kinds), {[0, 0], [1, 0], [0, 1], [1, 1]});
            neighbour = moment([1, 0], [1, 0]) - moment([1, 0], [1, 1]);
            expected = scale * k * (coincident(1) - coincident(2) - coincident(3) + 2 * coincident(4) + 2 * neighbour);

            testCase.verifyEqual(magneticXX(1, 1), expected, 'RelTol', 1e-10);
        end

        function testInfoMessageIsLogged(testCase)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.INFO, {sink});

            testCase.generate_(3, 3, 0.1, 0.1, 3);

            lines = sink.lines();
            testCase.verifyNumElements(lines, 1);
            testCase.verifySubstring(lines{1}, '[MomentAssembler.generate]: Assembled the generators of 12 rooftops');
        end

        function testPeakMemoryMatchesTheBudget(testCase)
            % Linux keeps the peak resident memory of the process and resets it on request.
            testCase.assumeTrue(isfile('/proc/self/clear_refs'), 'The peak memory is read from Linux /proc only.');

            % At the order 32 a lattice block of the far pass holds 1024 kernel values, 24 kB against the 560 bytes of
            % the small arrays of its offset, so their heap noise stays at a few percent (at the order 16 it is 9 %).
            % A budget of 256 MB holds 20 columns of the 513 blocks of a 512 x 512 grid.
            budget = 2^28;
            assemble = @(count) rooftopMoM.assembler.greensMatrix.MomentAssembler().generate( ...
                testCase.createContext(count, count, 0.05, 0.05, [32, 32], 2 * pi, budget));

            % A first assembly of the same size loads the code and the buffers of the matrix products, so only the
            % arrays of the generators and of the chunks remain to be measured.
            assemble(512);
            baseline = testCase.readProcessMemory('VmRSS');
            testCase.resetPeakMemory();
            generators = assemble(512);
            peak = testCase.readProcessMemory('VmHWM') - baseline;

            % rooftopMoM.Resources takes the generators and one chunk as the peak of the assembly. A far tile takes
            % whole columns of blocks, so it fills the chunk only to within one column, 5 % here. It counts its small
            % arrays in the chunk too (estimateBlockBytes_), since the heap of the C library may keep them resident
            % after they are freed; whether it does moves the peak by about 6 %.
            testCase.verifyLessThanOrEqual(peak, 1.02 * (budget + generators.memory()));
            testCase.verifyGreaterThanOrEqual(peak, 0.9 * budget);
        end
    end

    methods (Access = private)
        function generators = generate_(testCase, countX, countY, sizeX, sizeY, orders)
            context = testCase.createContext(countX, countY, sizeX, sizeY, orders, 2 * pi);
            generators = rooftopMoM.assembler.greensMatrix.MomentAssembler().generate(context);
        end
    end
end

function value = integrateAdaptively_(sizes, k, offset, kinds)
% integrateAdaptively_: Returns the moment M_alpha_beta, kinds = [alpha, beta], of exp(-ikR) / R at a touching offset,
%   integrated over the quadrants of [-1, 1]^2, which have the singular point at a corner.

    % W00, W10, W01 and W11 (Weights of the rooftop moments).
    weights = {@(nu) 1 - abs(nu), @(nu) (nu < 0) .* (1 + nu).^2 / 2 + (nu >= 0) .* (1 - nu.^2) / 2, ...
        @(nu) (nu < 0) .* (1 - nu.^2) / 2 + (nu >= 0) .* (1 - nu).^2 / 2, ...
        @(nu) (1 - abs(nu)).^2 .* (2 + abs(nu)) / 6};
    weight = weights{kinds(1) + 2 * kinds(2) + 1};
    distance = @(nu, eta) sqrt((sizes(1) * (nu - offset(1))).^2 + (sizes(2) * (eta - offset(2))).^2);
    integrand = @(nu, eta) weight(nu) .* (1 - abs(eta)) .* exp(-1i * k * distance(nu, eta)) ./ distance(nu, eta);
    ends = [-1, 0; 0, 1];
    value = 0;

    for m = 1:2
        for n = 1:2
            value = value + integral2(integrand, ends(m, 1), ends(m, 2), ends(n, 1), ends(n, 2), 'AbsTol', 0, ...
                'RelTol', 1e-13);
        end
    end
end
