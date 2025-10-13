function colormap = rainbow_desaturated_old(m)
% rainbow_desaturated(m) returns a desaturated rainbow colormap as an m x 3 matrix.
% m specifies the number of colors in the colormap.
% The colormap is desaturated compared to the standard MATLAB rainbow colormap.
if nargin < 1 || isempty(m)
   m = 256; % Default number of colormap entries
end
% Create the standard rainbow colormap
  standard_colormap = prism(m);

% Desaturate the colormap (example - can adjust this based on desired desaturation level)
  desaturated_colormap = 0.75 * standard_colormap; %  Reduce intensity to 75%

% Adjust red and green (optional) for better perception, particularly for colorblindness
  desaturated_colormap(1:m/2, 1) = desaturated_colormap(1:m/2, 1) * 0.7; %Reduce red
  desaturated_colormap(1:m/2, 2) = desaturated_colormap(1:m/2, 2) * 0.85; %Reduce green
  desaturated_colormap(m/2+1:m, 1) = desaturated_colormap(m/2+1:m, 1) * 0.85; %Reduce red
  desaturated_colormap(m/2+1:m, 2) = desaturated_colormap(m/2+1:m, 2) * 0.7; %Reduce green
%
  colormap = desaturated_colormap;
end