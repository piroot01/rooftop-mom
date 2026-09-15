classdef MemoryEfficientSplitFFTAcceleratedApplicator < rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator
% MemoryEfficientSplitFFTAcceleratedApplicator: Applies the Green's matrix using the split (subband) FFT with the
%   minimal peak memory. In 2D every contribution (subband, input component, output component) is transformed alone
%   in one grid-sized buffer and added to the output, so the peak is the output and a single grid array, the lower
%   bound of any split algorithm and a quarter of the working memory of the zero-padded FFT. It needs 64 FFTs
%   instead of the 24 of the SplitFFTAcceleratedApplicator. The single direction case is inherited, as the fastest
%   variant is already memory minimal there.

    properties (Constant, Access = private)
        % Number of column chunks copied between the grid buffer and the vectors, bounds the temporaries.
        chunkCount_ = 32;
    end

    methods
        function obj = MemoryEfficientSplitFFTAcceleratedApplicator(context, spectralMatrixData, rooftopContent)
            obj@rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator(context, spectralMatrixData, rooftopContent);
        end
    end

    methods
        function memory = workingMemory(obj)
        % workingMemory: Returns the peak memory of the arrays alive during one application (output included).
        %
        %   Output:
        %       memory: rooftopMoM.Memory of the peak working arrays.

            if obj.rooftopContent_ ~= rooftopMoM.basisFunctions.RooftopContent.BOTH
                memory = workingMemory@rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator(obj);
                return;
            end

            [fftSizeX, fftSizeY] = obj.spectralMatrixData_(1).blockSize();
            [rowCounts, columnCounts] = obj.componentLayout_();

            % The largest column chunk is copied up to three times at once (slice, output range and their sum).
            chunkElements = max(rowCounts .* ceil(columnCounts ./ min(obj.chunkCount_, columnCounts)));

            elements = obj.context_.rooftops().count() + fftSizeX * fftSizeY + 3 * chunkElements;
            memory = rooftopMoM.Memory(elements * ...
                rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.ComplexDouble).bytes(), 'B');
        end
    end

    methods (Access = protected)
        function output = applyBothDirs_(obj, input)
        % applyBothDirs_: Applies the Green's matrix of the x- and y-directed rooftops one contribution at a time.
        %
        %   Input:
        %       input: Nx1 vector of the rooftop coefficients, the x-directed rooftops first.
        %
        %   Output:
        %       output: Nx1 vector of the Green's matrix product.
        %
        %   NOTE: The copies between the grid buffer and the vectors are chunked column-wise, so no temporary of the
        %       grid size is ever created.

            [fftSizeX, fftSizeY] = obj.spectralMatrixData_(1).blockSize();
            [rowCounts, columnCounts, offsets] = obj.componentLayout_();

            columnEdges = cell(1, 2);
            for component = 1:2
                columnEdges{component} = round(linspace(0, columnCounts(component), ...
                    min(obj.chunkCount_, columnCounts(component)) + 1));
            end

            blockNames = {'blockXX', 'blockXY'; 'blockYX', 'blockYY'};

            output = zeros(numel(input), 1, 'like', 1i);
            buffer = zeros(fftSizeX, fftSizeY, 'like', 1i);

            for k = 1:4
                parityX = mod(k - 1, 2);
                parityY = (k - 1 - parityX) / 2;
                blocks = obj.spectralMatrixData_(k);

                for inputComponent = 1:2
                    rows = rowCounts(inputComponent);
                    edges = columnEdges{inputComponent};

                    for outputComponent = 1:2
                        % Copy the modulated input component into the buffer and zero its padding.
                        for chunk = 1:(numel(edges) - 1)
                            first = edges(chunk) + 1;
                            last = edges(chunk + 1);

                            part = reshape(input((offsets(inputComponent) + (first - 1) * rows + 1): ...
                                (offsets(inputComponent) + last * rows)), rows, []);

                            if parityX
                                part = part .* obj.modulationX_(1:rows);
                            end

                            if parityY
                                part = part .* obj.modulationY_(first:last);
                            end

                            buffer(1:rows, first:last) = part;
                        end
                        clear part;

                        buffer((rows + 1):end, :) = 0;
                        buffer(:, (columnCounts(inputComponent) + 1):end) = 0;

                        % Transform, mix and transform back in place.
                        buffer = fft(fft(buffer, [], 2), [], 1);
                        buffer = blocks.(blockNames{outputComponent, inputComponent})() .* buffer;
                        buffer = ifft(ifft(buffer, [], 1), [], 2);

                        % Add the demodulated output component to the output.
                        outputRows = rowCounts(outputComponent);
                        outputEdges = columnEdges{outputComponent};

                        for chunk = 1:(numel(outputEdges) - 1)
                            first = outputEdges(chunk) + 1;
                            last = outputEdges(chunk + 1);

                            part = buffer(1:outputRows, first:last);

                            if parityX
                                part = part .* obj.demodulationX_(1:outputRows);
                            end

                            if parityY
                                part = part .* obj.demodulationY_(first:last);
                            end

                            range = (offsets(outputComponent) + (first - 1) * outputRows + 1): ...
                                (offsets(outputComponent) + last * outputRows);
                            output(range) = output(range) + part(:);
                        end
                        clear part;
                    end
                end
            end
        end
    end

    methods (Access = private)
        function [rowCounts, columnCounts, offsets] = componentLayout_(obj)
        % componentLayout_: Returns the grid layout of the x- and y-directed rooftop components.
        %
        %   Outputs:
        %       rowCounts: 1x2 vector of the grid rows of the x- and y-directed components.
        %       columnCounts: 1x2 vector of the grid columns of the x- and y-directed components.
        %       offsets: 1x2 vector of the offsets of the components in the rooftop vectors.

            primitive = obj.context_.rooftops().primitive();
            countX = double(primitive.countX());
            countY = double(primitive.countY());

            rowCounts = [countX - 1, countX];
            columnCounts = [countY, countY - 1];
            offsets = [0, (countX - 1) * countY];
        end
    end
end
