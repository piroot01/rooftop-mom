function circulantKernel = generateRChanPreconditionerKernel(toeplitzKernel)
    if isempty(toeplitzKernel)
        circulantKernel = [];
        return;
    end

    [countX, countY] = size(toeplitzKernel);
    S = zeros(countX, countY, 'like', toeplitzKernel);
    S(2:end, 2:end) = rot90(toeplitzKernel(2:end, 2:end), 2);

    circulantKernel = toeplitzKernel + S;
end
