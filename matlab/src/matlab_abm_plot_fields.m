function matlab_abm_plot_fields(result, requested_times, fields, output_file)
%MATLAB_ABM_PLOT_FIELDS Plot saved fields using the papers' panel convention.

arguments
    result (1,1) struct
    requested_times (1,:) double
    fields (1,:) string
    output_file (1,1) string
end

indices = zeros(size(requested_times));
for index = 1:numel(requested_times)
    match = find(abs(result.time - requested_times(index)) < 1e-9, 1);
    if isempty(match)
        error('urbancrime:MissingSnapshot', ...
            'Time %g was not saved.', requested_times(index));
    end
    indices(index) = match;
end

single_time = isscalar(requested_times);
if single_time
    rows = 1;
    columns = numel(fields);
else
    rows = numel(fields);
    columns = numel(requested_times);
end
figure_handle = figure('Visible', 'off', 'Color', 'white');
layout = tiledlayout(rows, columns, 'TileSpacing', 'compact', 'Padding', 'compact');
for field_index = 1:numel(fields)
    field = fields(field_index);
    if ~isfield(result, field)
        error('urbancrime:UnknownField', 'Result has no field named %s.', field);
    end
    values = result.(field);
    limits = [min(values(:, :, indices), [], 'all'), ...
        max(values(:, :, indices), [], 'all')];
    if limits(1) == limits(2)
        limits(2) = limits(1) + eps(max(1, abs(limits(1))));
    end
    for time_index = 1:numel(indices)
        if single_time
            tile = field_index;
        else
            tile = (field_index - 1) * columns + time_index;
        end
        axis_handle = nexttile(layout, tile);
        imagesc(axis_handle, values(:, :, indices(time_index)));
        axis(axis_handle, 'image', 'xy', 'off');
        clim(axis_handle, limits);
        colormap(axis_handle, turbo(256));
        colorbar(axis_handle);
        if single_time
            title(axis_handle, field, 'Interpreter', 'none');
        elseif field_index == 1
            title(axis_handle, sprintf('t = %g', requested_times(time_index)));
        end
        if ~single_time && time_index == 1
            ylabel(axis_handle, field, 'Interpreter', 'none');
        end
    end
end
if single_time
    title(layout, sprintf('t = %g', requested_times(1)));
end
exportgraphics(figure_handle, output_file, 'Resolution', 220);
close(figure_handle);
end
