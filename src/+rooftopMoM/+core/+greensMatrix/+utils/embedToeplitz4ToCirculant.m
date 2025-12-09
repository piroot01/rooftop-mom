function H = embedToeplitz4ToCirculant(Tgen, cX, cY)
% embedToeplitz4ToCirculant: Circulant embedding of a 4th-order Toeplitz generator.
%
%   Purpose: Produce a reflected (mirrored) block suitable for circulant embedding
%       from a top-left Toeplitz generator block Tgen. The function accepts a
%       provided generator Tgen (which may be smaller or larger than the target
%       block) and returns a larger matrix H formed by reflecting rows and
%       columns of the prepared top-left block.
%
%   Inputs:
%       Tgen: numeric matrix
%           Provided generator block. Can be real or complex. May be any size;
%           it will be truncated or zero-padded to the target size (cY x cX).
%       cX, cY: positive integers
%           Target sizes for the top-left generator block in the Y (rows) and
%           X (columns) directions respectively.
%
%   Output:
%       H: numeric matrix
%           The circulant-embedded block constructed by reflecting the prepared
%           top-left block T (derived from Tgen). The intended shape depends on
%           cX and cY:
%               For the general case (cX>1 & cY>1):
%                   H = T( [1:cY, cY-1:-1:2], [1:cX, cX-1:-1:2] )
%               which yields a block of size (2*cY - 2) x (2*cX - 2).
%           Edge cases when one dimension equals 1 are handled explicitly
%           to produce "original then zeros" behaviour along that axis:
%               If cX == 1 and cY > 1 : H is (2*(cY-1)) x 2, with the first
%               column filled by reflected rows of T and the second column zeros.
%               If cY == 1 and cX > 1 : H is 2 x (2*(cX-1)), with the first
%               row filled by reflected columns of T and the second row zeros.
%
%   Additional notes: If Tgen is smaller than cY-by-cX it is zero-padded (bottom-right).
%       If Tgen is larger than cY-by-cX it is truncated to the top-left block.
%       The function attempts to preserve class and complexness where possible.
%       The central path uses zeros(...,'like',Tgen) so A inherits Tgen's class.
%       The special single-dimension branches create explicitly complex zeros
%       (complex(zeros(...))) to ensure complex-capable storage when needed.

    % Ignore empty kerneles.
    if isempty(Tgen)
        H = [];
        return;
    end

    % Sizes of provided Tgen (now interpret cX as rows and cY as columns)
    [cX0, cY0] = size(Tgen);

    % If cX and cY are empty set them same as cX0 and cY0
    if nargin < 2 || (isempty(cX) && isempty(cY))
        cX = cX0;
        cY = cY0;
    end

    % Prepare padded/truncated generator block T of size cX x cY.
    % Use 'like' to preserve numeric class (single/double, real/complex).
    if cX0 == cX && cY0 == cY
        % Tgen already matches target size: use directly.
        T = Tgen;
    else
        % Allocate minimal pad array and copy the top-left portion of Tgen.
        % This truncates if Tgen is larger, or zero-pads if Tgen is smaller.
        T = zeros(cX, cY, 'like', Tgen);

        cXMin = min(cX0, cX);
        cYMin = min(cY0, cY);

        % Copy the available entries from Tgen into T.
        T(1:cXMin, 1:cYMin) = Tgen(1:cXMin, 1:cYMin);
    end

    % ---------------------------------------------------------------------
    % Handle three distinct cases (roles swapped: cX = rows, cY = cols):
    %   1) cY == 1 & cX > 1   -> Produce (2*(cX-1)) x 2 block with first
    %                            column filled by reflected rows and second
    %                            column zeros.
    %   2) cY > 1 & cX == 1   -> Produce 2 x (2*(cY-1)) block with first row
    %                            filled by reflected columns and second row zeros.
    %   3) cX > 1 & cY > 1    -> General reflected 2D block using mirrored
    %                            index vectors for rows and columns.
    % ---------------------------------------------------------------------

    if cY == 1 && cX > 1
        % Column-single case (columns == 1): construct a tall, two-column matrix.
        % Use 'like' to preserve numeric class (single/double, real/complex).
        H = zeros(2 * (cX - 1), 2, 'like', Tgen);

        % Build reflected row index vector based on rows (cX):
        % [1,2,...,cX, cX-1,...,2]
        rowIndices = [1:cX, (cX - 1):(-1):2];

        % Fill the first column of H with the reflected column from T.
        % The second column remains zeros (original then zeros behaviour).
        H(:, 1) = T(rowIndices, :);
    elseif cY > 1 && cX == 1
        % Row-single case (rows == 1): construct a short, wide, two-row matrix.
        % Use 'like' to preserve numeric class (single/double, real/complex).
        H = zeros(2, 2 * (cY - 1), 'like', Tgen);

        % Build reflected column index vector based on columns (cY):
        % [1,2,...,cY, cY-1,...,2]
        columnIndices = [1:cY, (cY - 1):(-1):2];

        % Fill the first row of H with the reflected row from T.
        % The second row remains zeros (original then zeros behaviour).
        H(1, :) = T(:, columnIndices);
    elseif cX > 1 && cY > 1
        % General 2D reflection case: reflect both rows (cX) and columns (cY).
        % For rows: [1, 2, ..., cX, cX-1, ..., 2] (length 2*cX - 2)
        rowIndices = [1:cX, (cX - 1):(-1):2];

        % For cols: [1, 2, ..., cY, cY-1, ..., 2] (length 2*cY - 2)
        columnIndices = [1:cY, (cY - 1):(-1):2];

        % Select the reflected block directly using vectorized indexing.
        H = T(rowIndices, columnIndices);
    else
        % If both cX and cY are 1 the combination is unsupported by the
        % algorithm as written (no meaningful reflected block can be formed
        % under the current convention). Signal to the caller.
        error('Unsupported cX and cY parameters.');
    end
end
