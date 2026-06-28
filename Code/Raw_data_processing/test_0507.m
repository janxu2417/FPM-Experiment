close all;clear;clc;

%% ========== 1. 参数设置 ==========
% 低分辨率像素尺寸 (m) (对应你的 0.15NA 实际实验)
spsize  = 0.69e-6;
% % 上采样比例（假设 0.15NA 到 0.8NA 大约有 4~5 倍分辨率提升）
% pratio = 4;
% % 高分辨率像素尺寸 (m)
% psize = spsize / pratio;

% 低分辨率图像尺寸（像素）
m1 = 600;  % 高度
n1 = 600;  % 宽度

m = m1;
n = n1;

dark_noise_level = 0.005; % 固定本底噪声比例

% 光学参数
wlength = 6.28e-07; % 波长 (m)
NA  = 0.15;     % 模拟的低分辨率物镜数值孔径 (修改为你的实际物镜NA)
z   = 2e-6;        % 离焦量 (m) 模拟理想情况先设为0

%% ========== 2. 导入真实高分辨率图像并合成复振幅 ==========

num_img = 25;
im_high_HDR = zeros(m, n, num_img, 'double');
% im_med_HDR = zeros(m1, n1, num_img, 'double');

% 放在你的 MATLAB 脚本最前面
set(0, 'DefaultFigureVisible', 'on');

str_z = "+0.0430";
batch_folder = "0522";
raw_root = "Raw_input";

for i = 1:num_img
    % data_name = "0427\tif 原始格式\p" + i;
    % if (i <= 10)
    %     date_time = "0427\";
    % elseif (i <= 50)
    %     date_time = "0508\";
    % end
    % right_id = ceil(i / 10) * 10;
    % data_name = date_time + sprintf("%d-%d",right_id - 9, right_id) + "\图像_" + int2str(i);
    data_name = batch_folder + "\" + str_z + "\图像_" + int2str(i);
    data_dir = fullfile(raw_root, data_name + ".tif");
    % 读取你拍摄的 0.8NA TIFF 图像
    I_hr = imread(data_dir);
    I_hr = double(I_hr(:,:,1));
    fprintf('图 %s 的尺寸为 %d %d\n', data_name, size(I_hr));
    fprintf('这幅图的物理最大光强值是: %d\n', max(I_hr(:)));
    % figure;
    % imshow(I_hr, []);
    % title('MATLAB 自动对比度显示的真实数据');

    % 为了演示代码能跑通，如果找不到图，生成一个测试矩阵（你用时将其注释掉）
    %I_hr = double(imresize(rgb2gray(imread('cameraman.tif')), [4096, 3000]));

    % 0. 归一化 + 扣除相机固定本底噪声
    I_hr = I_hr ./ max(I_hr(:));
    I_hr = I_hr - dark_noise_level;
    I_hr(I_hr < 0) = 0;

    % 2. 截取中心区域 (m x n) 对应 FPM 的视场
    center_y = round(size(I_hr, 1) / 2);
    center_x = round(size(I_hr, 2) / 2);
    row_start = center_y + 1;
    row_end = center_y + m;
    col_start = center_x - n/2 + 1;
    col_end = center_x + n/2;
    I_crop = I_hr(row_start:row_end, col_start:col_end);

    % 3. 从强度计算振幅 (Amplitude = sqrt(Intensity))
    amp = sqrt(I_crop);

    % 4. 物理模拟：合成"真实"的相位
    % 在生物细胞等透明样品中，相位延迟通常正比于样本厚度或吸收（振幅）
    % 这里我们假设一个随振幅变化的弱相位物体（最大相位差 pi）
    % phase = (1 - amp) * pi;
    phase = 0;
    % % 合成高分辨率基准复振幅 O
    O = amp .* exp(1j * phase);

    im_high_HDR(:,:,i) = amp;
    % im_med_HDR(:,:,i) = imresize(amp, 0.5, "bicubic", "Antialiasing", true);

    % figure('Name', '基准真值 Ground Truth (Ideal)');
    % imshow(abs(O), []); title('HR 振幅 (真实数据)');
    % subplot(122); imshow(im_med_HDR(:,:,i), []); title('HR 振幅 (压缩一半)');

    %     % subplot(122); imshow(angle(O), []); title('HR 相位 (合成物理模型)');
    % folder = fullfile('D:','1_徐前','物理','傅里叶叠层显微术','Code', 'Raw_input', '0508', '高分辨原图');
    % filename = sprintf('LED_%d_高分辨率原图.png', i);
    % saveas(gcf, fullfile(folder, filename));

    % himFT = fftshift(fft2(I_crop));
    % figure('Name', '高分辨率频谱分析');
    % % 使用 log(1 + abs(...)) 进行动态范围压缩
    % % subplot(121);
    % imshow(log(1 + abs(himFT(50:end-50, 50:end-50))), []);
    % title('Log Amplitude (对数振幅谱)');

    % subplot(122);
    % imshow(angle(himFT(50:end-50, 50:end-50)), []);
    % title('Phase (相位谱)');
    % folder = fullfile('D:','1_徐前','物理','傅里叶叠层显微术','Code', 'Raw_input', '0427');
    % folder = fullfile('D:','1_徐前','物理','傅里叶叠层显微术','Code', 'Raw_input', '0508', '高分辨频谱');
    % filename = sprintf('LED_%d_高分辨率频谱.png', i);
    % saveas(gcf, fullfile(folder, filename));
end
out_name = "FPM_RawData_z " + str_z + "_2.mat";
folder = fullfile('D:','1_徐前','物理','傅里叶叠层显微术','Code', 'Raw_input', batch_folder);
save(fullfile(folder, out_name), 'z', 'NA', 'wlength', 'spsize', "im_high_HDR");


%{
%% Prepare the experimental data
% Add necessary folders into the current working directory
% Load data file
data_name = 'USAF_red';
data_dir = ['Data\' data_name '.mat'];
load(data_dir); % refer to 'data_description.txt' for more details
% Display raw images
is_show = 'all'; % 'center' shows the first low-res raw image; 'all' dynamically shows all low-res images
if strcmp(is_show,'center')
    figure(1);
    set(gcf,'outerposition',get(0,'ScreenSize'))
    imshow(imlow_HDR(:,:,1),[]);
    title(['raw image ' num2str(1)]);
elseif strcmp(is_show,'all')
    for k = 1:size(imlow_HDR,3)
        figure(1);
        set(gcf,'outerposition',get(0,'ScreenSize'))
        imshow(imlow_HDR(:,:,k),[]);
        title(['raw image ' num2str(k)]); pause(0.1);
    end
end

%% Set up the experiment parameters
xstart = 18; ystart = 20; % absolute coordinate of initial LED
arraysize = 3;% side length of lit LED array
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);
H      = 40.24; % distance between LEDs and sample, in mm
LEDp   = 4;     % distance between adjacent LEDs, in mm
nglass = 1.52;  % refraction index of glass substrate
t      = 1;     % glass thickness, in mm
[kx, ky, NAt] = k_vector(xlocation-xstart, ylocation-ystart, H, LEDp, nglass, t, theta, xint, yint, arraysize^2);

%% Reconstruct by FP algorithm
NA          = 0.1;      % objective NA
spsize      = 1.845e-6; % pixel size of low-res image on sample plane, in m
upsmp_ratio = 4;        % upsampling ratio
psize       = spsize/upsmp_ratio; % pixel size of high-res image on sample plane, in m

opts.loopnum    = 10;   % iteration number
opts.alpha      = 1;    % '1' for ePIE, other value for rPIE
opts.beta       = 1;    % '1' for ePIE, other value for rPIE
opts.gamma_obj  = 1;    % the step size for object updating
opts.gamma_p    = 1;    % the step size for pupil updating
opts.eta_obj    = 0.2;  % the step size for adding momentum to object updating
opts.eta_p      = 0.2;  % the step size for adding momentum to pupil updating
opts.T          = 1;    % do momentum every T images. '0' for no momentum during the recovery; integer, generally (0, arraysize^2].
opts.aberration = aberration; % pre-calibrated aberration, if available

used_idx = 1:1:arraysize^2; % choose which raw image is used, for example, 1:2:arraysize^2 means do FPM recovery with No1 image, No3 image, No5 image......
imlow_used = imlow_HDR(1:1488,301:1788,used_idx);
kx_used = kx(used_idx);
ky_used = ky(used_idx);
[him, tt, fprobe, imlow_HDR1] = himrecover(imlow_used, kx_used, ky_used, NA, wlength, spsize, psize, z, opts);
% 检查 him 中的异常值
nan_count = sum(isnan(him(:)));
inf_count = sum(isinf(him(:)));
zero_count = sum(him(:) == 0);

fprintf('NaN 数量: %d\n', nan_count);
fprintf('Inf 数量: %d\n', inf_count);
fprintf('零值数量: %d\n', zero_count);

figure;
set(gcf,'outerposition',get(0,'ScreenSize'))
subplot(121);imshow(abs(him(50:end-50,50:end-50)),[]);title('Amplitude');
subplot(122);imshow(angle(him(50:end-50,50:end-50)),[]);title('Phase');
disp(['Wavelength: ',num2str(wlength.*1e+9),' nm, Loop: ',num2str(opts.loopnum)]);
disp(['Maximum illumination NA = ',num2str(max(NAt(used_idx)))]);

%% Save the results
out_dir = 'Results';
mkdir(out_dir); addpath(out_dir);
out_name = [data_name '_result.mat'];
save([out_dir,'\',out_name],'him','fprobe','tt','imlow_HDR1');
%}