function [X, Y] = generate_grid(x_min, y_min, num_rows, num_cols, ell)
%

%

x      = x_min + ((1 : num_rows) - 1 + 0.5) * ell;
y      = y_min + ((1 : num_cols) - 1 + 0.5) * ell;
[X, Y] = meshgrid(x, y);
end