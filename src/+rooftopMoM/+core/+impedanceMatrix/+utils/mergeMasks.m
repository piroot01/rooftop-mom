function [composed, shifts] = mergeMasks(maskXX, maskXY, maskYY)
% mergeMasks: Combines sub-masks into a composite block mask.
%   It assembles a composite mask by arranging input masks into a 2x2 block format.
%   The off-diagonal 'maskXY' is shifted to avoid index collision.
%
%   Inputs:
%       maskXX: Numerical index mask for the top-left (XX) block.
%       maskXY: Numerical index mask for the off-diagonal (XY) block.
%       maskYY: Numerical index mask for the bottom-right (YY) block.
%
%   Outputs:
%       composed: The resulting composite 2x2 block mask.
%       shifts: A struct containing the offset value applied to 'maskXY'.

    % Determine the required shift from the maximum index in the XX mask.
    maxXX = max(maskXX(:));
    maskYYShifted = maskYY;

    % Shift the XY mask to ensure its indices sequentially follow the XX indices.
    maskXYShifted = rooftopMoM.core.impedanceMatrix.utils.shiftMask(maskXY, maxXX);

    % Assemble the final block matrix using the original and shifted masks.
    composed = [maskXX, maskXYShifted; ...
                maskXYShifted.', maskYYShifted];

    % Record the shift amount used for the XY block for external reference.
    shifts.xy = maxXX;
end
