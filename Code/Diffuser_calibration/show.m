%% ===== Read tif from a Raw_input folder, display and save as png =====
clear; clc; close all;

% 1. Input folder
input_folder = fullfile( ...
    'D:\1_徐前\物理\傅里叶叠层显微术\Code\Raw_input', ...
    '0529\光强修正');
    % '0605', ...
    % '毛玻璃01');

output_folder = fullfile(input_folder, 'png_export');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 2. Find tif files and sort by numeric suffix
files = dir(fullfile(input_folder, '图像_*.tif'));
if isempty(files)
    error('No tif files found in: %s', input_folder);
end

file_nums = zeros(numel(files), 1);
for i = 1:numel(files)
    token = regexp(files(i).name, '图像_(\d+)\.tif', 'tokens', 'once');
    if isempty(token)
        error('Unexpected filename: %s', files(i).name);
    end
    file_nums(i) = str2double(token{1});
end

[~, order] = sort(file_nums);
files = files(order);
file_nums = file_nums(order);

% 3. Export each tif to png
for i = 1:numel(files)
    tif_path = fullfile(files(i).folder, files(i).name);
    I = imread(tif_path);

    if ndims(I) >= 3
        I_show = I(:, :, 1);
    else
        I_show = I;
    end

    fig = figure('Color', 'w', 'Visible', 'on');
    imshow(I_show, []);
    title(sprintf('File: %s', files(i).name), 'Interpreter', 'none');

    png_name = sprintf('图像_%d.png', file_nums(i));
    exportgraphics(gca, fullfile(output_folder, png_name), 'Resolution', 150);

    close(fig);
end

fprintf('Saved %d png files to:\n%s\n', numel(files), output_folder);

% 4. Optional: create a montage overview
enable_montage_overview = false;

if enable_montage_overview
    num_show = min(numel(files), 25);
    img_stack = cell(1, num_show);

    for i = 1:num_show
        tif_path = fullfile(files(i).folder, files(i).name);
        I = imread(tif_path);
        if ndims(I) >= 3
            I = I(:, :, 1);
        end
        img_stack{i} = mat2gray(I);
    end

    figure('Color', 'w', 'Position', [100, 100, 1200, 900]);
    montage(img_stack, 'Size', [ceil(sqrt(num_show)), ceil(sqrt(num_show))]);
    title(sprintf('Montage of first %d tif images', num_show), 'Interpreter', 'none');

    frame = getframe(gcf);
    imwrite(frame.cdata, fullfile(output_folder, 'montage_first_images.png'));
end