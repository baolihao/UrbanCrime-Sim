function smooth_cmap = rainbow_desaturated()
% RGB control points for ParaView's Rainbow Desaturated
cmap = [
    0.2081, 0.1663, 0.5292;
    0.0060, 0.3522, 0.5524;
    0.0381, 0.5041, 0.4316;
    0.0986, 0.6250, 0.2650;
    0.2906, 0.7109, 0.0199;
    0.5172, 0.7684, 0.0000;
    0.7411, 0.7969, 0.0000;
    0.9020, 0.8000, 0.0000;
    1.0000, 0.8000, 0.0000;
    1.0000, 0.7250, 0.0000;
    1.0000, 0.6460, 0.0000;
];

% Interpolate to create smoother colormap with 256 entries
npts = 256;
x = linspace(0, 1, size(cmap,1));
xq = linspace(0, 1, npts);
smooth_cmap = interp1(x, cmap, xq, 'pchip');
end