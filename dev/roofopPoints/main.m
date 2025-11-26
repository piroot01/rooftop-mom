clc; clearvars;

profile on;

rng(0);

countX = 3;
countY = 3;
order = 2;
repeatCount = 1;

orderSquared = order^2;

singleQuadraturePointsPositive = randi(10, [3, orderSquared]);

tic;
shifts = generateRectGridShifts(countX, countY, 1, 1);
toc;

count = size(shifts, 2);
cols = orderSquared * count;

for n = 1:repeatCount
    tic;
    shiftsDistributed = repelem(shifts, 1, orderSquared);
    quadraturePoints1 = repmat(singleQuadraturePointsPositive, 1, count) + shiftsDistributed
    toc;
    clear quadraturePoints1;

    tic;
    qP = reshape(singleQuadraturePointsPositive, 3, orderSquared, 1);
    s = reshape(shifts, 3, 1, count);
    quadraturePoints2 = reshape(s + qP, 3, cols)
    toc;
    clear quadraturePoints2;
end

profile viewer;

function shifts = generateRectGridShifts(countX, countY, sizeX, sizeY, idx)
% generateRectGridShifts: Generates shift vectors for a rectangular grid.
%
%   shifts = generateRectGridShifts(..., idx)
%       -> returns shifts only for the linear indices in idx. idx may be:
%            - a vector of indices (e.g. 4:16)
%            - a two-element vector [start end] (treated as start:end)
%            - a single scalar
%
% Notes:
%   - Linear indices follow the same ordering as reshape(ndgrid(...),1,[]): 
%     the x-index (countX) varies fastest, y-index (countY) slower.
%   - Returned shifts are doubles: 3 x numel(idx).
%
% Example:
%   % shifts for columns 4..16 only
%   s = generateRectGridShifts(2000,1000,1.0,1.0, 4:16);

    % total number of cells
    N = double(countX) * double(countY);

    % if no idx provided, keep original behaviour (small grids ok)
    if nargin < 5 || isempty(idx)
        [xn, yn] = ndgrid(double(0:(countX - 1)) * double(sizeX), ...
                          double(0:(countY - 1)) * double(sizeY));
        shifts = [reshape(xn, 1, []); reshape(yn, 1, []); zeros(1, numel(xn))];
        return
    end

    % normalize idx:
    if numel(idx) == 2 && idx(1) < idx(2) && all(mod(idx,1)==0) % treat as range [start end]
        idx = idx(1):idx(2);
    end
    idx = idx(:).';  % row vector

    % convert to 0-based linear offsets
    k0 = double(idx) - 1;

    % compute 0-based grid coordinates (vectorized, no big temporary grids)
    ix = mod(k0, double(countX));        % 0..countX-1 (x index)
    iy = floor(k0 / double(countX));     % 0..countY-1 (y index)

    % scale by cell size and build shifts
    xcoords = ix * double(sizeX);
    ycoords = iy * double(sizeY);

    shifts = [xcoords; ycoords; zeros(1, numel(idx))];
end