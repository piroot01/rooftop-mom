classdef (Abstract) Base < matlab.unittest.TestCase
% Base: Common fixture of the unit tests.
%   It puts src on the path, so a test class also runs with only the folder tests on the path. The shared logger and
%   the resources are reset before and after every test, so no test depends on the state left by another one.

    methods (TestClassSetup)
        function addSourcePath(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
        end
    end

    methods (TestMethodSetup)
        function resetLogger(testCase)
            logger = rooftopMoM.logging.Logger.instance();
            logger.reset();
            testCase.addTeardown(@() logger.reset());
        end

        function resetResources(testCase)
            rooftopMoM.Resources.reset();
            testCase.addTeardown(@() rooftopMoM.Resources.reset());
        end
    end

    methods (Access = protected)
        function bytes = readProcessMemory(~, field)
        % readProcessMemory: Returns the resident memory of the process, 'VmRSS', or its peak, 'VmHWM', from the Linux
        %   file /proc/self/status.

            text = fileread('/proc/self/status');
            value = regexp(text, [field, ':\s*(\d+) kB'], 'tokens', 'once');
            bytes = 1024 * str2double(value{1});
        end

        function resetPeakMemory(~)
        % resetPeakMemory: Sets the peak resident memory of the process back to the current one; Linux does so when 5
        %   is written to /proc/self/clear_refs.

            fileId = fopen('/proc/self/clear_refs', 'w');
            fprintf(fileId, '5');
            fclose(fileId);
        end

        function rooftops = createRooftops(~, countX, countY, sizeX, sizeY, order)
        % createRooftops: Returns the rooftops of a grid with the quadrature of the given order.

            primitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
            quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
            rooftops = rooftopMoM.basisFunctions.UniformRooftops(primitive, quadrature);
        end

        function context = createContext(testCase, countX, countY, sizeX, sizeY, orders, wavenumber, memoryBytes)
        % createContext: Returns the AssemblerContext of a grid with the quadrature orders [near, far], or one order
        %   for both, by default with a 64 MB chunk budget. The rooftops get the far order, as in Problem.

            if nargin < 8
                memoryBytes = 2^26;
            end

            orders = [orders(1), orders(end)];
            rooftops = testCase.createRooftops(countX, countY, sizeX, sizeY, orders(2));
            context = rooftopMoM.assembler.AssemblerContext(rooftops, wavenumber, memoryBytes, orders);
        end

        function product = createProduct(~, name, rooftops, generators)
        % createProduct: Returns the product with the Green's matrix of the applicator class name (without the
        %   package) built from the spectra it needs.

            primitive = rooftops.primitive();

            if strcmp(name, 'FFTApplicator')
                spectra = rooftopMoM.assembler.greensMatrix.utils.generateEmbeddingSpectra(generators, ...
                    primitive.sizeX(), primitive.sizeY());
            else
                spectra = rooftopMoM.assembler.greensMatrix.utils.generateSubbandSpectra(generators, ...
                    primitive.sizeX(), primitive.sizeY());
            end

            product = feval(['rooftopMoM.assembler.greensMatrix.', name], rooftops, spectra, generators.content());
        end

        function matrix = generateDenseMatrix(~, generators, rooftops)
        % generateDenseMatrix: Returns the dense Green's matrix built from its generators by the definition: the block
        %   Toeplitz vector potential part minus D.' * E * D, with E the face interaction matrix and D the divergence
        %   of the rooftops (1 / faceSize on the first face of a rooftop, -1 / faceSize on the second).

            primitive = rooftops.primitive();
            gridSize = [primitive.countX(), primitive.countY()];
            [faceI, faceJ] = ndgrid(0:(gridSize(1) - 1), 0:(gridSize(2) - 1));
            faceMatrix = toeplitz_(generators.electric(), [faceI(:).'; faceJ(:).']);

            directions = 'xy';
            magneticParts = {generators.magneticXX(), generators.magneticYY()};
            faceSizes = [primitive.sizeX(), primitive.sizeY()];
            magneticBlocks = cell(1, 2);
            divergence = zeros(prod(gridSize), 0);

            for d = 1:2
                count = rooftops.countDir(directions(d));
                first = rooftops.faceOffsets(directions(d), 1:count);
                second = first + rooftops.referenceRooftop(directions(d)).secondFaceOffset;

                magneticBlocks{d} = toeplitz_(magneticParts{d}, first);

                faces = [sub2ind(gridSize, first(1, :) + 1, first(2, :) + 1), ...
                    sub2ind(gridSize, second(1, :) + 1, second(2, :) + 1)];
                values = [ones(1, count), -ones(1, count)] / faceSizes(d);
                divergence = [divergence, full(sparse(faces, [1:count, 1:count], values, prod(gridSize), count))]; %#ok<AGROW>
            end

            matrix = blkdiag(magneticBlocks{:}) - divergence.' * faceMatrix * divergence;
        end
    end
end

function matrix = toeplitz_(generator, offsets)
    % The generators are even in both offsets, so an entry depends on the absolute offset difference only.
    if isempty(offsets)
        matrix = zeros(0);
        return;
    end

    matrix = generator(sub2ind(size(generator), abs(offsets(1, :).' - offsets(1, :)) + 1, ...
        abs(offsets(2, :).' - offsets(2, :)) + 1));
end
