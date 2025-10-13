function highway_val = get_highway(X, Y, ell)
% function highway_val = get_highway(X, Y, ell)

% the highway function in Eqn (4.1) and (4.2) on page 13 in the note

highway_val      = exp(-20 * (X + Y - ell).^2);
ind              = Y >= ell/2;
highway_val(ind) = 0;
end