%% 模拟FP成像过程：生成高分辨率原图 + 多角度低分辨率图像 + FP_recover_code所需输入
clear; close all; clc;

%% ========== 1. 参数设置 ==========
% 低分辨率像素尺寸 (m)
spsize  = 1.845e-6;   
% 高分辨率像素尺寸 (m)
psize = 1.845e-6 / 4;  % 假设高分辨率像素很小（超分辨）
% 上采样比例（整数）
pratio = round(spsize / psize);  

% 低分辨率图像尺寸（像素）
m1 = 201;  % 高度
n1 = 201;  % 宽度
% 高分辨率图像尺寸
m = m1 * pratio;
n = n1 * pratio;

% 光学参数
wlength = 6.28e-07; % 波长 (m)
NA_obj  = 0.1;  % 物镜数值孔径
z       = 2.5e-5;    % 离焦量 (m)
%{
%% ========== 2. 生成随机高分辨率复振幅图像 ==========
 振幅：随机物体（0~1）
amp = rand(m, n);
amp = amp ./ max(amp(:));  % 归一化
 相位：随机相位（0~2π）
phase = 2 * 0.0 * rand(m, n);
 复振幅
O = amp .* exp(1j * phase);

intensityMatrix = imresize(rgb2gray(imread('Human Skin.jpg')), [m, n]);
amp = double(intensityMatrix) .^ 0.5;
amp = amp ./ max(amp(:));
%}
%{
addpath(genpath(pwd));
data_name = 'MouseKidney_green_result';
data_dir = ['Results\' data_name '.mat']; 
load(data_dir);

% 复振幅
O = him;

% 保存高分辨率原图
figure;
subplot(121); imshow(abs(O(50:end-50,50:end-50)), []); title('HR 振幅');
subplot(122); imshow(angle(O(50:end-50,50:end-50)), []); title('HR 相位');
%save('HR_object.mat', 'O', 'psize');
%}
%% ========== 3. 计算每个LED的照明波矢 ==========
xstart = 18; ystart = 20; % absolute coordinate of initial LED
arraysize = 15; % side length of lit LED array
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);
H      = 90.88; % distance between LEDs and sample, in mm
LEDp   = 4;     % distance between adjacent LEDs, in mm
nglass = 1.52;  % refraction index of glass substrate
t      = 1;     % glass thickness, in mm
xint   = 0.0;
yint   = 0.0;
theta  = 0.0;
[kx_norm, ky_norm, NAt] = k_vector(xlocation-xstart, ylocation-ystart, H, LEDp, nglass, t, theta, xint, yint, arraysize^2);
% 调用 k_vector.m 计算归一化波矢 (kx, ky)
LED_num = arraysize^2;
% 转换为实际波矢 (1/μm)
k0 = 2 * pi / wlength;
kx = kx_norm * k0;
ky = ky_norm * k0;
dkx = 2*pi/(psize*n); dky = 2*pi/(psize*m); % 波矢最小间隔

%pratio = round(spsize/psize); % upsampling ratio 上采样比例（低分辨率/高分辨率像素尺寸比）
%m = pratio*m1; n = pratio*n1;
%NAfilx = NA*(1/wlength)*n*psize; NAfily = NA*(1/wlength)*m*psize; % m1*spsize = m*psize 截止频率半径
%kmax = pi/psize; % the max wave vector of the OTF
%kx2 = -kmax:kmax/((n-1)/2):kmax; ky2 = -kmax:kmax/((m-1)/2):kmax; % odd N x,y方向波矢采样点
%[kxm, kym] = meshgrid(kx2,ky2); % 生成二维波矢网格

%% ========== 4. 生成瞳孔函数 P(kx, ky) ==========
[KX, KY] = meshgrid( (-n/2 : n/2 - 1) * (2*pi/(n*psize)), ...
                     (-m/2 : m/2 - 1) * (2*pi/(m*psize)) );
K = sqrt(KX.^2 + KY.^2);
P = (K <= NA_obj * k0);  % 圆形瞳孔

%% ========== 5. 多角度成像模拟 ==========
%imseqlow = zeros(m1, n1, LED_num);  % 存储低分辨率图像序列

imseqlow1 = zeros(m1, n1, LED_num);

for i = 1:LED_num
    % 平面波照明
    [x, y] = meshgrid( (1:n) * psize, (1:m) * psize );
    illu = exp(1j * (kx(1,i) * x + ky(1,i) * y));
    kxc=round((n+1)/2-kx(1,1)/dkx);
    kyc=round((m+1)/2-ky(1,1)/dky);
    kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
    kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
    % 物面场分布
    %U = O .* illu .* tt(1,i+(15-1)*LED_num);
    %U = O .* illu;
    %O = O ./ tt(1,i+9*LED_num);
    % 傅里叶变换
    U_ft = fftshift(fft2(O .* illu));

    %U_ft1 = fftshift(fft2(O));
    
    % 瞳孔滤波

    %V_ft=U_ft.*P;

    O_j=U_ft(kyl:kyh,kxl:kxh);
    % 提取当前瞳孔对应的傅里叶域数据（O_j = 高分辨率傅里叶变换的局部区域）
    V_ft1=O_j.*fprobe; % low-pass filter

    % 逆傅里叶变换到像面
    %V = ifft2(ifftshift(V_ft));

    V1 = ifft2(ifftshift(V_ft1));
    % 取强度 + 下采样得到低分辨率图像
    %I_low  = imresize(abs(V), 1/pratio, 'bilinear');
    
    I_low1 = abs(V1);
    % 存储
    %imseqlow(:, :, i) = I_low;

    imseqlow1(:, :, i) = I_low1;
end

%% ========== 6. 保存 FP_recover_code.m 所需输入数据 ==========
save('FP_input_data.mat', ...
     'imseqlow1', ...  % 低分辨率图像序列 (m1 x n1 x LED_num)
     'kx_norm', 'ky_norm', ...  % 每个LED的波矢 (1 x LED_num)
     'spsize', ...     % 低分辨率像素尺寸 (μm)
     'psize', ...    % 高分辨率像素尺寸 (μm)
     'wlength', ...   % 波长 (μm)
     'NA_obj', ...    % 物镜NA
     'z');            % 离焦量 (μm)

fprintf('数据生成完成！保存文件：FP_input_data.mat 和 HR_object.mat\n');

opts.loopnum    = 10;   % iteration number
opts.alpha      = 1;    % '1' for ePIE, other value for rPIE
opts.beta       = 1;    % '1' for ePIE, other value for rPIE
opts.gamma_obj  = 1;    % the step size for object updating
opts.gamma_p    = 1;    % the step size for pupil updating
opts.eta_obj    = 0.2;  % the step size for adding momentum to object updating
opts.eta_p      = 0.2;  % the step size for adding momentum to pupil updating
opts.T          = 1;    % do momentum every T images. '0' for no momentum during the recovery; integer, generally (0, arraysize^2].
opts.aberration = 0;    % pre-calibrated aberration, if available

is_show = 'YES';
if strcmp(is_show,'YES')
    for k = 1:100:201
        figure();
        set(gcf,'outerposition',get(0,'ScreenSize'));
        %set(gcf,'outerposition',get(0,'ScreenSize'))
        subplot(121); imshow(imlow_HDR1(:,:,k), []); title(['raw image given' num2str(k)]);
        %imshow(20*log10(abs(fftshift(fft2(imlow_HDR1(:,:,k)))))); title(['raw image given' num2str(k)]);
        %subplot(121); imshow(30*log10(abs(fftshift(fft2(imseqlow(:,:,k)))))); title(['raw image generated P ' num2str(k)]);
        %subplot(122); imshow(10*log10(abs(fftshift(fft2(imseqlow1(:,:,k)))))); title(['raw image generated f ' num2str(k)]);
        
        %imshow(imseqlow(:,:,k),[]); title(['raw image generated using P ' num2str(k)]);
        subplot(122); imshow(imseqlow1(:,:,k),[]); title(['raw image generated using fprobe ' num2str(k)]);
    end
end
FP_test_recover(imseqlow1, kx_norm, ky_norm, NAt, psize, spsize, wlength, NA_obj, z, opts, O);