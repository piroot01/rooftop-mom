function points = generateQuadraturePointsSlice(singlePoints, shifts)
    singlePoints = reshape(singlePoints, 3, [], 1);
    shifts = reshape(shifts, 3, 1, []);
    points = reshape(shifts + singlePoints, 3, []);
end
