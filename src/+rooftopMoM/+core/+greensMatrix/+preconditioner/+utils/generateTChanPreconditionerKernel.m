function circulantKernel = generateTChanPreconditionerKernel(toeplitzKernel)
    if isempty(toeplitzKernel)
        circulantKernel = [];
        return;
    end

    [countX, countY] = size(toeplitzKernel);
    circulantKernel = zeros(countX, countY, 'like', toeplitzKernel);
    count = countX * countY;

    for p = 0:(countX - 1)
        ip1 = p + 1;
        ip2 = mod(countX - p, countX) + 1;

        for q = 0:(countY - 1)
            jq1 = q + 1;
            jq2 = mod(countY - q, countY) + 1;

            % weights
            w00 = (countX - p) * (countY - q);
            w10 = p * (countY - q);
            w01 = (countX - p) * q;
            w11 = p * q;

            circulantKernel(ip1, jq1) = (w00 * toeplitzKernel(ip1, jq1) ...
                + w10 * toeplitzKernel(ip2, jq1) ...
                + w01 * toeplitzKernel(ip1, jq2) ...
                + w11 * toeplitzKernel(ip2, jq2)) / count;
        end
    end
end
