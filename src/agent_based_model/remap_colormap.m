function [cMap, theMap] = remap_colormap(A, kind)
%

if nargin < 2, kind = 1; end
switch kind
  case 1
    theMap  = colormap('jet');
  case 2
    theMap  = rainbow_desaturated;
  otherwise
end
A_max       = max(max(A));
A_min       = min(min(A));
[m, n]      = size(A);
if A_max ~= A_min
  A_scaled  = (A - A_min)/(A_max - A_min);
else
  A_scaled  = ones(m, n);
end
nRows_max   = 246;
nRows_min   = 6;
color_ind   = round((nRows_max - nRows_min) * A_scaled + nRows_min);
color_ind   = color_ind(:);
cMap        = theMap(color_ind, :);
cMap        = reshape(cMap, [m, n, 3]);
end