function generators = generate(obj, context)
% generate: Returns the generators of the rooftop directions present on the grid.
%
%   Input:
%       context: rooftopMoM.assembler.AssemblerContext.
%
%   Output:
%       generators: rooftopMoM.assembler.greensMatrix.Generators.

    start = tic;
    rooftops = context.rooftops();
    primitive = rooftops.primitive();
    counts = [primitive.countX(), primitive.countY()];
    sizes = [primitive.sizeX(), primitive.sizeY()];
    wavenumber = context.wavenumber();
    orders = context.quadratureOrders();

    % The moments are divided by the squared cell area; scales turns them into entries of the electric and of the
    % magnetic generators.
    scales = 1i * rooftopMoM.assembler.Constants.Z0 / (4 * pi) * prod(sizes)^2 * [1 / wavenumber, wavenumber];
    parts = allocateParts_(rooftops);

    % The near pairs have the offsets 0 to nearCounts - 1 in each direction.
    nearCounts = countNearOffsets_(counts, sizes, obj.nearGap_);

    % The second argument is the number of parts of [-1, 0] and [0, 1] in x and y, not a size: the rounded aspect
    % ratio, e.g. [3, 1] for a 3 x 1 cell. The near rule thus splits the longer cell side into parts of about the
    % shorter one, so the near singularity of an elongated pair lies as many parts away from the nodes as for square
    % cells. The far pairs are smooth and need no split.
    nearRule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(orders(1), round(sizes / min(sizes)));
    farRule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(orders(2), [1, 1]);

    logger = rooftopMoM.logging.Logger.instance();
    logger.debug('Cell offsets: %d near with the quadrature order %d, %d far with the order %d.', ...
        prod(nearCounts), orders(1), prod(counts) - prod(nearCounts), orders(2));

    parts = integrateCellPairs_(parts, context, scales, nearRule, nearCounts, [0, 0]);
    parts = integrateCellPairs_(parts, context, scales, farRule, counts, nearCounts);

    % The touching pairs integrated the remainder of the extraction; the extracted part follows in closed form.
    touching = [0, 1, 0, 1; 0, 0, 1, 1];
    inside = all(touching < counts(:), 1);
    extracted = rooftopMoM.assembler.singularIntegrals.generateExtractedMoments(sizes(1), sizes(2), wavenumber);
    parts = rooftopMoM.assembler.greensMatrix.utils.accumulateMoments(parts, extracted(:, inside), ...
        touching(:, inside), scales);

    generators = rooftopMoM.assembler.greensMatrix.Generators(parts.magneticXX, parts.magneticYY, parts.electric, ...
        rooftops.content());

    logger.info('Assembled the generators of %d rooftops in %.3f s (%.1f MB).', rooftops.count(), toc(start), ...
        generators.memory() / 2^20);
end
