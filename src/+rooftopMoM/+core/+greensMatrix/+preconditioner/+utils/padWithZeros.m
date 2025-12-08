function output = padWithZeros(input, cX, cY)
    if isempty(input)
        output = [];
        return;
    end

    [r0, c0] = size(input);
    if r0 == cX && c0 == cY
        output = input;
        return;
    end

    output = zeros(cX, cY, 'like', input);

    r = 1:min(r0, cX);
    c = 1:min(c0, cY);

    output(r, c) = input(r, c);
end
