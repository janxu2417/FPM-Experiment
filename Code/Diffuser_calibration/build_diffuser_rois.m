function roi_info = build_diffuser_rois(image_size, main_roi_size, aux_roi_size, aux_roi_offsets, center_xy)
%BUILD_DIFFUSER_ROIS Build one central ROI and several auxiliary ROIs.
% image_size: [height, width]
% main_roi_size / aux_roi_size: scalar or [height, width]
% aux_roi_offsets: [num_aux x 2], each row is [dx, dy] in pixels
% center_xy: optional [x, y] center in pixels. If omitted, image center is used.

image_h = image_size(1);
image_w = image_size(2);
if nargin < 5 || isempty(center_xy)
    center_x = round((image_w + 1) / 2);
    center_y = round((image_h + 1) / 2);
else
    if ~(isnumeric(center_xy) && numel(center_xy) == 2)
        error('center_xy must be empty or a two-element numeric vector [x, y].');
    end
    center_xy = double(center_xy(:)).';
    center_x = round(center_xy(1));
    center_y = round(center_xy(2));
end

main_size = local_parse_roi_size(main_roi_size);
aux_size = local_parse_roi_size(aux_roi_size);

main_rect = local_build_rect(center_x, center_y, main_size, image_w, image_h, 'main ROI');
roi_info = struct();
roi_info.image_size = [image_h, image_w];
roi_info.center_xy = [center_x, center_y];
roi_info.main = local_build_roi_struct('main', [0, 0], main_rect);

num_aux = size(aux_roi_offsets, 1);
aux = repmat(local_build_roi_struct('', [0, 0], [0, 0, 0, 0]), num_aux, 1);
for k = 1:num_aux
    dx = aux_roi_offsets(k, 1);
    dy = aux_roi_offsets(k, 2);
    rect = local_build_rect(center_x + dx, center_y + dy, aux_size, image_w, image_h, ...
        sprintf('aux ROI %d', k));
    aux(k) = local_build_roi_struct(sprintf('aux_%d', k), [dx, dy], rect);
end
roi_info.aux = aux;
roi_info.all = [roi_info.main; roi_info.aux(:)];

end

function roi_size = local_parse_roi_size(value)
if isscalar(value)
    roi_size = [value, value];
elseif isnumeric(value) && numel(value) == 2
    roi_size = reshape(value, 1, 2);
else
    error('ROI size must be a scalar or a two-element vector.');
end

roi_size = round(roi_size);
if any(roi_size <= 0)
    error('ROI size must be positive.');
end
end

function rect = local_build_rect(center_x, center_y, roi_size, image_w, image_h, roi_name)
half_w = floor((roi_size(2) - 1) / 2);
half_h = floor((roi_size(1) - 1) / 2);
x1 = center_x - half_w;
x2 = x1 + roi_size(2) - 1;
y1 = center_y - half_h;
y2 = y1 + roi_size(1) - 1;

if x1 < 1 || x2 > image_w || y1 < 1 || y2 > image_h
    error('%s is out of image bounds.', roi_name);
end

rect = [x1, x2, y1, y2];
end

function roi = local_build_roi_struct(name, offset_xy, rect)
roi = struct();
roi.name = name;
roi.offset_xy = offset_xy;
roi.rect = rect;
roi.width = rect(2) - rect(1) + 1;
roi.height = rect(4) - rect(3) + 1;
end
