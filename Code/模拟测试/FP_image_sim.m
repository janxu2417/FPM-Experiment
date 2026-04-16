%% 模拟FP成像过程：使用真实的0.8NA高分辨率原始TIFF图像
clear; close all; clc;

%% ========== 1. 参数设置 ==========
% 低分辨率像素尺寸 (m) (对应你的 0.15NA 实际实验)
spsize  = 1.845e-6;   
% 上采样比例（假设 0.15NA 到 0.8NA 大约有 4~5 倍分辨率提升）
pratio = 5;  
% 高分辨率像素尺寸 (m)
psize = spsize / pratio;  

% 低分辨率图像尺寸（像素）
m1 = 400;  % 高度
n1 = 400;  % 宽度
% 需要的高分辨率图像尺寸
m = m1 * pratio;
n = n1 * pratio;

% 光学参数
wlength = 6.28e-07; % 波长 (m)
NA_obj  = 0.15;     % 模拟的低分辨率物镜数值孔径 (修改为你的实际物镜NA)
z       = 0;        % 离焦量 (m) 模拟理想情况先设为0

%% ========== 2. 导入真实高分辨率图像并合成复振幅 ==========
% 读取你拍摄的 0.8NA TIFF 图像
I_hr = imread('Raw_input\图像_high.tif'); 
% 【注】：这里我用伪代码代替，请换成你实际的文件名
I_hr = double(I_hr(:,:,1));

fprintf('这幅图的物理最大光强值是: %d\n', max(I_hr(:))); 
figure;
imshow(I_hr, []); 
title('MATLAB 自动对比度显示的真实数据');
% 为了演示代码能跑通，如果找不到图，生成一个测试矩阵（你用时将其注释掉）
%I_hr = double(imresize(rgb2gray(imread('cameraman.tif')), [4096, 3000])); 

% 0. 扣除相机固定本底噪声
I_hr = I_hr - 3;
I_hr(I_hr < 0) = 0;
% 1. 归一化强度到 [0, 1]
I_hr = I_hr ./ max(I_hr(:));

% 2. 截取中心区域 (m x n) 对应 FPM 的视场
center_y = round(size(I_hr, 1) / 2);
center_x = round(size(I_hr, 2) / 2);
I_crop = I_hr(center_y - m/2 + 1 : center_y + m/2, center_x -n/4 + 1 : center_x + n*3/4);

% 3. 从强度计算振幅 (Amplitude = sqrt(Intensity))
amp = sqrt(I_crop);

% 4. 物理模拟：合成"真实"的相位
% 在生物细胞等透明样品中，相位延迟通常正比于样本厚度或吸收（振幅）
% 这里我们假设一个随振幅变化的弱相位物体（最大相位差 pi）
phase = (1 - amp) * pi;  

% 合成高分辨率基准复振幅 O
O = amp .* exp(1j * phase);

figure('Name', '基准真值 Ground Truth (0.8 NA)');
subplot(121); imshow(abs(O), []); title('HR 振幅 (真实数据)');
subplot(122); imshow(angle(O), []); title('HR 相位 (合成物理模型)');

%% ========== 3. 计算每个LED的照明波矢 ==========
xstart = 18; ystart = 20; % absolute coordinate of initial LED
arraysize = 10; % side length of lit LED array
% 【注】LED_location 和 k_vector 是你自己的函数，保持不变
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);
H      = 40.24; % mm
LEDp   = 4;     % mm
nglass = 1.52;  
t      = 1;     
xint   = 0.0; yint = 0.0; theta  = 1.0;
[kx_norm, ky_norm, NAt] = k_vector(xlocation-xstart, ylocation-ystart, H, LEDp, nglass, t, theta, xint, yint, arraysize^2);

LED_num = arraysize^2;
k0 = 2 * pi / wlength;
kx = kx_norm * k0;
ky = ky_norm * k0;
dkx = 2*pi/(psize*n); 
dky = 2*pi/(psize*m);

%% ========== 4. 生成 0.15 NA 的瞳孔函数 P (低分辨率尺寸) ==========
% 注意：在频域裁剪法中，瞳孔的大小应该是 m1 x n1 (低分辨率的像素数)
[KX, KY] = meshgrid( (-n1/2 : n1/2 - 1) * (2*pi/(n1*spsize)), ...
                     (-m1/2 : m1/2 - 1) * (2*pi/(m1*spsize)) );
K = sqrt(KX.^2 + KY.^2);
% 理想圆形低通滤波器
P = double(K <= NA_obj * k0);  

%% ========== 5. 严谨的物理多角度成像模拟 (正向模型) ==========
imseqlow1 = zeros(m1, n1, LED_num);

% 提前将高分辨率基准原图进行傅里叶变换
O_ft = fftshift(fft2(O));

fprintf('开始生成 0.15 NA 模拟采集序列...\n');
for i = 1:LED_num
    % 1. 计算当前 LED 对应的频域偏移中心坐标
    kxc = round((n+1)/2 - kx(1,i)/dkx);
    kyc = round((m+1)/2 - ky(1,i)/dky);
    
    % 2. 确定频域裁剪的边界 (大小为 m1 x n1)
    kyl = round(kyc - (m1-1)/2);  kyh = round(kyc + (m1-1)/2);
    kxl = round(kxc - (n1-1)/2);  kxh = round(kxc + (n1-1)/2);
    
    % 防止高频大角度越界报错的保护机制
    if kyl < 1 || kxl < 1 || kyh > m || kxh > n
        warning('LED %d 角度过大，超出高分网格边界！全置0', i);
        imseqlow1(:, :, i) = 0;
        continue;
    end
    
    % 3. 提取当前孔径截取的频域数据
    O_j = O_ft(kyl:kyh, kxl:kxh);
    
    % 4. 乘以低通瞳孔函数（物理上的物镜低通滤波）
    V_ft = O_j .* P; 
    
    % 5. 逆傅里叶变换到相机靶面
    V = ifft2(ifftshift(V_ft));
    
    % 6. 相机采集：取振幅 (因为你的 himrecover 默认输入是振幅，若输入强度需平方)
    imseqlow1(:, :, i) = abs(V); 
end
fprintf('序列生成完成！\n');

%% ========== 6. 保存与后续恢复接口 ==========
% 保存数据供其他程序调用
save('FP_input_data_Sim.mat', 'imseqlow1', 'kx_norm', 'ky_norm', 'spsize', 'psize', 'wlength', 'NA_obj', 'z');

opts.loopnum    = 10;   % iteration number
opts.alpha      = 1;    % '1' for ePIE, other value for rPIE
opts.beta       = 1;    % '1' for ePIE, other value for rPIE
opts.gamma_obj  = 1;    % the step size for object updating
opts.gamma_p    = 1;    % the step size for pupil updating
opts.eta_obj    = 0.2;  % the step size for adding momentum to object updating
opts.eta_p      = 0.2;  % the step size for adding momentum to pupil updating
opts.T          = 1;    % do momentum every T images. '0' for no momentum during the recovery; integer, generally (0, arraysize^2].
opts.aberration = 0;    % pre-calibrated aberration, if available

% 这里你可以直接调用你的恢复函数（代码略）
%FP_test_recover(imseqlow1, kx_norm, ky_norm, NAt, spsize, wlength, NA_obj, z, opts, O);

% 可视化检查：对比 0.8 NA 的基准图和 0.15 NA 正入射拍到的图
figure('Name', '系统分辨率对比');
subplot(121); imshow(abs(O), []); title('基准 0.8 NA (HR)');
subplot(122); imshow(imseqlow1(:, :, ceil(LED_num/2)), []); title('模拟中心正入射 0.15 NA (LR)');