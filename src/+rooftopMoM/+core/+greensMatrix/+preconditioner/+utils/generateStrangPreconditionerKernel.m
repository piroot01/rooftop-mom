function circulantKernel = generateStrangPreconditionerKernel(toeplitzKernel)
    if isempty(toeplitzKernel)
        circulantKernel = [];
        return;
    end

    [countX, countY] = size(toeplitzKernel);
    circulantKernel = zeros(countX, countY, 'like', toeplitzKernel);

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
