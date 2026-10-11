function circulant = foldGenerator(generator, window, even)
% foldGenerator: Returns the generator of a circulant from a Toeplitz generator: every offset is weighted by a window
%   and added to its offset modulo the period. The circulant approximations of the kernels are such folds and differ
%   in the window only.
%
%   Inputs:
%       generator: (2 * countX - 1) x (2 * countY - 1) two-sided generator with the offset 0 in the middle, or with
%           even the countX x countY one-sided generator of the offsets 0, ..., count - 1.
%       window: Function handle window(offsets, count, dimension) that returns the weights of the offsets
%           -(count - 1), ..., count - 1 along the dimension 1 (x) or 2 (y).
%       even: Optional true if the generator is one-sided and even, g(-p, q) = g(p, -q) = g(p, q) (default false).
%
%   Output:
%       circulant: countX x countY generator of the circulant.
%
%   NOTE: In one direction, a 4 x 4 Toeplitz matrix T with the offsets t(-3), ..., t(3) gets the circulant C of the
%       period 4, whose first column c(0), ..., c(3) is the fold:
%           T = [t(0) t(-1) t(-2) t(-3)        C = [c(0) c(3) c(2) c(1)
%                t(1) t(0)  t(-1) t(-2)             c(1) c(0) c(3) c(2)
%                t(2) t(1)  t(0)  t(-1)             c(2) c(1) c(0) c(3)
%                t(3) t(2)  t(1)  t(0)]             c(3) c(2) c(1) c(0)]
%       c(0) = w(0) * t(0) and c(k) = w(k) * t(k) + w(k - 4) * t(k - 4) for k = 1, 2, 3. The window of T. Chan,
%       w(p) = (4 - |p|) / 4, gives c(1) = (3 * t(1) + t(-3)) / 4: the offset 1 occurs three times in T, its partner
%       -3 once. In two directions the fold runs along x and then along y. An even generator has t(k - 4) = t(4 - k),
%       so its one-sided generator t(0), ..., t(3) is enough.

    if nargin > 2 && even
        [countX, countY] = size(generator);
    else
        countX = (size(generator, 1) + 1) / 2;
        countY = (size(generator, 2) + 1) / 2;
        even = false;
    end

    weightsX = reshape(window(-(countX - 1):(countX - 1), countX, 1), [], 1);
    weightsY = reshape(window(-(countY - 1):(countY - 1), countY, 2), 1, []);

    if even
        % The offset k keeps its weight w(k) and gains its partner k - count, which the even generator holds at
        % count - k, with the weight w(k - count). The offset 0 has no partner, so row 1 stands in with the weight 0.
        partnersX = [0; weightsX(1:(countX - 1))];
        partnersY = [0, weightsY(1:(countY - 1))];
        circulant = weightsX(countX:end) .* generator + partnersX .* generator([1, countX:(-1):2], :);
        circulant = circulant .* weightsY(countY:end) + circulant(:, [1, countY:(-1):2]) .* partnersY;
    else
        % The window is the product of a weight along x and a weight along y.
        weighted = weightsX .* generator .* weightsY;

        % Fold along x, the rows; the transpose turns the columns into rows for the fold along y.
        circulant = foldRows_(weighted);
        circulant = foldRows_(circulant.').';
    end
end

function folded = foldRows_(generator)
    % The offset p sits in the row p + count. The offsets 0, ..., count - 1 keep their place, the offsets
    % -(count - 1), ..., -1 are added to the offsets 1, ..., count - 1.
    count = (size(generator, 1) + 1) / 2;
    folded = generator(count:end, :);
    folded(2:end, :) = folded(2:end, :) + generator(1:(count - 1), :);
end
