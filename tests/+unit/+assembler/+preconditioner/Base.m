classdef (Abstract) Base < unit.Base
% Base: Common fixture of the tests of the circulant preconditioners, dense two-level Toeplitz matrices and the
%   Fourier matrix of a grid, with the x index running fastest.

    methods (Access = protected)
        function matrix = generateToeplitzMatrix(~, generator)
        % generateToeplitzMatrix: Returns the dense two-level Toeplitz matrix of a two-sided generator.

            grid = (size(generator) + 1) / 2;
            [i, j] = ndgrid(0:(grid(1) - 1), 0:(grid(2) - 1));
            matrix = generator(sub2ind(size(generator), i(:) - i(:).' + grid(1), j(:) - j(:).' + grid(2)));
        end

        function fourier = generateFourierMatrix(~, grid)
        % generateFourierMatrix: Returns the unitary two-dimensional Fourier matrix of a grid.

            fourier = kron(fft(eye(grid(2))), fft(eye(grid(1)))) / sqrt(prod(grid));
        end
    end
end
