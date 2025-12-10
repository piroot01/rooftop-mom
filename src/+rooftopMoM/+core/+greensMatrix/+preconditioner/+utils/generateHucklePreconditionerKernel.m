function circulantKernel = generateHucklePreconditionerKernel(toeplitzKernel, pX, pY, centerScale)
    if isempty(toeplitzKernel)
        circulantKernel = [];
        return;
    end

    if nargin < 4
        centerScale = 1;
    end

    [countX, countY] = size(toeplitzKernel);
    circulantKernel = zeros(countX, countY, 'like', toeplitzKernel);

    % safe iteration limits (matches intent of original, but prevents OOB)
    maxI = min(pX - 1, countX - 1);
    maxJ = min(pY - 1, countY - 1);

    % small local aliases for speed
    T = toeplitzKernel;
    C = circulantKernel;
    invPX = 1.0 / pX;
    invPY = 1.0 / pY;

    % center
    C(1,1) = centerScale * T(1,1);

    % first row contributions (i = 0)
    for j = 1:maxJ
        jp = j + 1;
        jneg = countY - j + 1;
        wY = 1.0 - j * invPY;
        v = T(1, jp) * wY;
        C(1, jp) = C(1, jp) + v;
        C(1, jneg) = C(1, jneg) + v;
    end

    % first column and interior blocks (i > 0)
    for i = 1:maxI
        ip = i + 1;
        ineg = countX - i + 1;
        wX = 1.0 - i * invPX;

        % first-column contributions (j = 0)
        vcol = T(ip, 1) * wX;
        C(ip, 1) = C(ip, 1) + vcol;
        C(ineg, 1) = C(ineg, 1) + vcol;

        % interior (j > 0)
        for j = 1:maxJ
            jp = j + 1;
            jneg = countY - j + 1;
            wY = 1.0 - j * invPY;
            w = wX * wY;
            v = T(ip, jp) * w;
            C(ip, jp) = C(ip, jp) + v;
            C(ip, jneg) = C(ip, jneg) + v;
            C(ineg, jp) = C(ineg, jp) + v;
            C(ineg, jneg) = C(ineg, jneg) + v;
        end
    end

    % return result
    circulantKernel = C;
end

%{
function circulantKernel = generateHucklePreconditionerKernel(toeplitzKernel, pX, pY)
    if isempty(toeplitzKernel)
        circulantKernel = [];
        return;
    end

    [countX, countY] = size(toeplitzKernel);
    circulantKernel = zeros(countX, countY, 'like', toeplitzKernel);

    idxX = (1:(pX - 1)).';
    idxY = (1:(pY - 1));

    wX = 1 - (idxX ./ pX);
    wY = 1 - (idxY ./ pY);

    idxXPos = idxX + 1;
    idxXNeg = countX - idxX + 1;

    idxYPos = idxY + 1;
    idxYNeg = countY - idxY + 1;

    toeplitzKernelBlock = toeplitzKernel(idxX + 1, idxY + 1);
    values = (wX .* wY) .* toeplitzKernelBlock;

    % center
    circulantKernel(1, 1) = toeplitzKernel(1, 1);

    % first column (j=0): use toeplitzKernel(idxXPos,1) for both +i and -i
    if ~isempty(idxX)
        colvals = toeplitzKernel(idxXPos, 1) .* wX;
        circulantKernel(idxXPos, 1) = circulantKernel(idxXPos, 1) + colvals;
        circulantKernel(idxXNeg, 1) = circulantKernel(idxXNeg, 1) + colvals;
    end

    % first row (i=0): use toeplitzKernel(1, idxYPos) for both +j and -j
    if ~isempty(idxY)
        rowvals = toeplitzKernel(1, idxYPos) .* wY;
        circulantKernel(1, idxYPos) = circulantKernel(1, idxYPos) + rowvals;
        circulantKernel(1, idxYNeg) = circulantKernel(1, idxYNeg) + rowvals;
    end

    % four quadrants for i>0, j>0: add same values to all symmetric positions
    if ~isempty(values)
        circulantKernel(idxXPos, idxYPos) = circulantKernel(idxXPos, idxYPos) + values;
        circulantKernel(idxXPos, idxYNeg) = circulantKernel(idxXPos, idxYNeg) + values;
        circulantKernel(idxXNeg, idxYPos) = circulantKernel(idxXNeg, idxYPos) + values;
        circulantKernel(idxXNeg, idxYNeg) = circulantKernel(idxXNeg, idxYNeg) + values;
    end

    circulantKernel
end
%}
