function twoSided = extendEvenly(generator)
% extendEvenly: Returns the two-sided generator of a Toeplitz matrix that is even in each offset,
%   g(-p, q) = g(p, -q) = g(p, q) and so also g(-p, -q) = g(p, q), from its one-sided generator of p, q >= 0. The
%   generators of the Green's function are such, as they depend on |p| and |q| only.
%
%   Input:
%       generator: countX x countY generator of the offsets 0, ..., count - 1 (may be empty).
%
%   Output:
%       twoSided: (2 * countX - 1) x (2 * countY - 1) generator with the offset 0 in the middle.

    [countX, countY] = size(generator);
    twoSided = generator([countX:(-1):2, 1:countX], [countY:(-1):2, 1:countY]);
end
