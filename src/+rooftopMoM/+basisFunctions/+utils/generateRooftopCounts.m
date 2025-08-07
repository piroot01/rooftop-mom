function [c1, c2] = generateRooftopCounts(n1, n2)
% generateRooftopCounts: Compute rooftop counts for a rectangular grid.
%   It returns the number of rooftop basis functions oriented along the two principal axes of a fully-populated rectangular mesh.
%
%   Inputs:
%       n1: Number of nodes along the 1st axis.
%       n2: Number of nodes along the 2nd axis.
%
%   Outputs:
%       c1: Number of rooftop elements spanning x-directed edges.
%       c2: Mumber of rooftop elements spanning y-directed edges.

    c1 = (n1 - 1) * n2;
    c2 = (n2 - 1) * n1;
end

