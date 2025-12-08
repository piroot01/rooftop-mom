function circulantKernel = generateStrangPreconditionerKernel(toeplitzKernel, centerScale)
    if isempty(toeplitzKernel)
        circulantKernel = [];
        return;
    end

    if nargin < 2
        centerScale = 1;
    end

    [countX, countY] = size(toeplitzKernel);
    circulantKernel = zeros(countX, countY, 'like', toeplitzKernel);

    % Adjust the most dominant kernel value.
    toeplitzKernel(1, 1) = toeplitzKernel(1, 1) * centerScale;

    halfX = floor(countX / 2);
    halfY = floor(countY / 2);

    for nX = 0:(countX - 1)
        if nX <= halfX
            p = nX + 1;
        else
            p = countX - nX + 1;
        end

        for nY = 0:(countY - 1)
            if nY <= halfY
                q = nY + 1;
            else
                q = countY - nY + 1;
            end

            circulantKernel(nX + 1, nY + 1) = toeplitzKernel(p, q);
        end
    end
end
