function bytes = estimateBlockBytes_(order)
% estimateBlockBytes_: Returns the peak bytes per lattice block of integrateLattice_, a block of order x order kernel
%   values: 24 bytes per value while the distances and the kernel exist together (a real and a complex double), plus
%   560 bytes for the small arrays that form and add the moments of the block's offset.
%   Those small arrays hold up to 35 complex values per offset at once: the moments, the sums of a quadrant and the
%   parts added in integrateTile_, then the moments, their combinations and the terms in accumulateMoments. They live
%   after the kernel is freed, but the heap of the C library keeps arrays below 32 MiB resident once they are freed,
%   so in the next tile they sit next to the kernel. Measured: 445 bytes per block at the order 16. Next to the 24 *
%   order^2 bytes of the kernel they are 7 % at the order 16 and 27 % at the order 8, too much to leave out.

    bytes = 24 * order^2 + 560;
end
