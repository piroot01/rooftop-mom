function output = applyBothDirections_(obj, input)
% applyBothDirections_: Applies the Green's matrix of the x- and y-directed rooftops one contribution at a time.
%
%   Input:
%       input: Nx1 vector of the rooftop coefficients, the x-directed rooftops first.
%
%   Output:
%       output: Nx1 vector of the product.
%
%   NOTE: The copies between the grid buffer and the vectors are chunked column-wise, so no temporary of the grid
%       size is ever created.

    [fftSizeX, fftSizeY] = obj.spectra_(1).gridSize();
    [rowCounts, columnCounts, offsets] = obj.componentLayout_(obj.rooftops_);

    columnEdges = cell(1, 2);
    for component = 1:2
        columnEdges{component} = round(linspace(0, columnCounts(component), ...
            min(obj.chunkCount_, columnCounts(component)) + 1));
    end

    blockNames = {'blockXX', 'blockXY'; 'blockYX', 'blockYY'};

    output = zeros(numel(input), 1, 'like', 1i);
    buffer = zeros(fftSizeX, fftSizeY, 'like', 1i);

    for k = 1:4
        % Subband k has the odd frequencies in x for k = 2, 4 and in y for k = 3, 4.
        isOddX = mod(k - 1, 2) == 1;
        isOddY = k > 2;
        blocks = obj.spectra_(k);

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

                    if isOddX
                        part = part .* obj.modulationX_(1:rows);
                    end

                    if isOddY
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

                    if isOddX
                        part = part .* obj.demodulationX_(1:outputRows);
                    end

                    if isOddY
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
