function [him,tt,fmaskpro,imseqlow]=himrecover(imseqlow,kx,ky,NA,wlength,spsize,psize,z,opts)
% FP algorithm to recover high-resolution image from low-resolution measured images
% Input:
%       imseqlow: low-res measurements, [m1 x n1 x numim] matrix
%       kx,ky: normalized wavevector of each LED, which should times k0 later 归一化波矢数组
%       NA: objective NA
%       wlength: central peak wavelength of LED, in m
%       spsize: pixel size of low-res image on sample plane, in m
%       psize: pixel size of high-res image on sample plane, in m
%       z: known defocus distance, in m
%       opts: other optimization parameters
% Output:
%       him: recovered high-res image, [m x n] matrix
%       tt: recorded intensity ratio between estimated and measured low-res amplitude, used for intensity correction
%       fmaskpro: recovered pupil function
%       imseqlow: low-res amplitudes after intensity correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% input arguments check
if ~isfield(opts,'loopnum')
    opts.loopnum = 10;
end
if ~isfield(opts,'alpha')
    opts.alpha = 1;
end
if ~isfield(opts,'beta')
    opts.beta = 1;
end
if ~isfield(opts,'gamma_obj')
    opts.gamma_obj = 1;
end
if ~isfield(opts,'gamma_p')
    opts.gamma_p = 1;
end
if ~isfield(opts,'eta_obj')
    opts.eta_obj = 0;
end
if ~isfield(opts,'eta_p')
    opts.eta_p = 0;
end
if ~isfield(opts,'T')
    opts.T = 0;
end
if ~isfield(opts,'aberration')
    opts.aberration = 0;
end
loopnum = opts.loopnum;
alpha = opts.alpha;
beta = opts.beta;
gamma_obj = opts.gamma_obj;
gamma_p = opts.gamma_p;
eta_obj = opts.eta_obj;
eta_p = opts.eta_p;
T = opts.T;
aberration = opts.aberration;

%% k-space parameterization
[m1, n1, numim] = size(imseqlow);
pratio = round(spsize/psize); % upsampling ratio 上采样比例（低分辨率/高分辨率像素尺寸比）
m = pratio*m1; n = pratio*n1;
k0 = 2*pi/wlength;
kx = k0*kx; ky = k0*ky;
NAfilx = NA*(1/wlength)*n*psize; NAfily = NA*(1/wlength)*m*psize; % m1*spsize = m*psize 截止频率半径
kmax = pi/psize; % the max wave vector of the OTF
dkx = 2*pi/(psize*n); dky = 2*pi/(psize*m); % 波矢最小间隔
kx2 = -kmax:kmax/((n-1)/2):kmax; ky2 = -kmax:kmax/((m-1)/2):kmax; % odd N x,y方向波矢采样点
[kxm, kym] = meshgrid(kx2,ky2); % 生成二维波矢网格
kzm = sqrt(k0^2-kxm.^2-kym.^2);

%% prior knowledge of aberration
H2 = exp(1j.*z.*real(kzm)).*exp(-abs(z).*abs(imag(kzm))); % define the defocus aberration if it is known or you want to test it
astigx = 0; astigy = 0; % define the astigmatism aberration if it is known or you want to test it
[M1, N1] = meshgrid(1:m1,1:n1);
zn = astigx*gzn(max(m1,n1),2*max(round(NAfily),round(NAfilx)),2,2)+...
     astigy*gzn(max(m1,n1),2*max(round(NAfily),round(NAfilx)),-2,2);
zn = imresize(zn,[m1,n1]);
if  any(aberration ~= 0,'all')
    fmaskpro = aberration; % pre-calibrated aberrations 若有预校准像差，直接使用
else
    fmaskpro = 1.*double(((M1-(m1+1)/2)/NAfily).^2+((N1-(n1+1)/2)/NAfilx).^2<=1)... % low-pass filter
    .*H2(round((m+1)/2-(m1-1)/2):round((m+1)/2+(m1-1)/2),round((n+1)/2-(n1-1)/2):round((n+1)/2+(n1-1)/2))... % defocus aberration
    .*exp(pi*1j.*zn); % astigmatism aberration 像散像差（zn为Zernike函数生成的像散项）
    % In this example, we can test the effect of the defocus aberration (z) and astigmatism aberrations (astigx, astigy) 
    % If the aberration is unknown, one can test different z, astigx, astigy for the best result
    % Gradient descent can also be used to update z, astigx, astigy (did not implement here)
    % Higher-order Zernike modes can be tested in a similar manner
end
%% initialization
him = imresize(sum(imseqlow,3),[m,n]); 
himFT = fftshift(fft2(him));

%% main part to optimize estimate of high-res image
for i = 1:2 % 2轮初始迭代 initial iterations to get a rough estimate of high-res image
    for i3 = 1:numim
        if max(max(imseqlow(:,:,i3))) == 0
            continue
        end
        % when the image size is even, there will be a half pixel displacement for the center. 
        kxc=round((n+1)/2-kx(1,i3)/dkx);  
        kyc=round((m+1)/2-ky(1,i3)/dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
        %fprintf('%d %d %d %d\n',kxl, kxh, m, n);
        O_j=himFT(kyl:kyh,kxl:kxh);
        % 提取当前瞳孔对应的傅里叶域数据（O_j = 高分辨率傅里叶变换的局部区域）
        lowFT=O_j.*fmaskpro; % low-pass filter
        im_lowFT=ifft2(ifftshift(lowFT)); % 逆傅里叶变换到空间域，得到估计的低分辨率图像

        updatetemp=pratio^2.*imseqlow(:,:,i3); % 采样比例校正（强度守恒）
        im_lowFT=updatetemp.*exp(1j.*angle(im_lowFT)); 
        lowFT_p=fftshift(fft2(im_lowFT)); % 反馈更新高分辨率傅里叶域数据（对应交替投影步骤）

        himFT(kyl:kyh,kxl:kxh)=himFT(kyl:kyh,kxl:kxh)+...
            conj(fmaskpro)./(max(max((abs(fmaskpro)).^2))).*(lowFT_p-lowFT);
        % conj 共轭; abs().^2 透过率振幅平方; max(max()) 二维矩阵所有元素的最大值
    end
end
countimg = 0;
tt = ones(1,loopnum*numim);

% for momentum method
vobj0 = zeros(m,n);
vp0 = zeros(m1, n1);
ObjT = himFT; 
PT = fmaskpro;

for i = 1:loopnum % 迭代次数 
    for i3 = 1:numim 
        countimg=countimg+1; 
        % 改进 1：全零图/极弱信号图跳过逻辑 (防止高频擦除与瞳孔污染)
        % =================================================================
        if max(max(imseqlow(:,:,i3))) == 0  % 或者设定一个极小的值如 < 1e-4
            % 该角度无有效物理信号，直接跳过本轮更新
            % 不更新频谱、不更新瞳孔、不更新光强校正系数
            
            % 如果启用了动量法，仍需保证 countimg 达到 T 时进行动量更新，
            % 但不影响当前跳过的逻辑
            if countimg == T
                vobj = eta_obj.*vobj0 + (himFT - ObjT);  himFT = ObjT + vobj;
                vobj0 = vobj;  ObjT = himFT;
                vp = eta_p.*vp0 + (fmaskpro - PT);  fmaskpro = PT + vp;
                vp0 = vp;  PT = fmaskpro;
                countimg = 0;
            end
            continue; % 直接进入下一个 LED 角度
        end
        % =================================================================
        kxc=round((n+1)/2-kx(1,i3)/dkx);  
        kyc=round((m+1)/2-ky(1,i3)/dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
        O_j=himFT(kyl:kyh,kxl:kxh); 
        lowFT=O_j.*fmaskpro;
        im_lowFT=ifft2(ifftshift(lowFT));
        % LED intensity correction 
        mean_imseqlow = mean(mean(pratio^2*abs(imseqlow(:,:,i3))));
        if mean_imseqlow > 0.05 % 设定一个基于你图像归一化范围的安全下限
            tt(1,i3+(i-1)*numim)=(mean(mean(abs(im_lowFT))))/mean_imseqlow; 
        else
            fprintf('i3 = %d, mean_imseqlow = %d \n', i3, mean_imseqlow);
        end
        if i>2
            imseqlow(:,:,i3)=imseqlow(:,:,i3).*tt(1,i3+(i-1)*numim); 
        end     
        % 反馈更新高分辨率傅里叶域数据
        updatetemp=pratio^2.*imseqlow(:,:,i3);
        im_lowFT=updatetemp.*exp(1j.*angle(im_lowFT)); 
        lowFT_p=fftshift(fft2(im_lowFT));
        % 1. 目标图像更新：用误差信号优化高分辨率频谱，alpha为平衡系数，避免分母过小
        % 检查更新量
        himFT(kyl:kyh,kxl:kxh) = himFT(kyl:kyh,kxl:kxh) +...
            gamma_obj.*conj(fmaskpro).*(lowFT_p-lowFT)./...
                ((1-alpha).*abs(fmaskpro).^2 + alpha.*max(max(abs(fmaskpro).^2)));

        % 2. 瞳孔函数更新：同步优化瞳孔像差，beta为平衡系数，对应"联合样本-瞳孔恢复"
        % 优化后的瞳孔更新逻辑
        illum_NA = sqrt(kx(1,i3)^2 + ky(1,i3)^2) / k0; % 当前 LED 的照明 NA
        if illum_NA < NA * 0.9 % 仅用明场及明暗场交界处的数据来算像差
            fmaskpro = fmaskpro + gamma_p.*conj(O_j).*((lowFT_p-lowFT))./...
                ((1-beta).*abs(O_j).^2 + beta.*max(max(abs(O_j).^2)));   
        end

        if mod(i3, 20) == 0
            fprintf('迭代状态: i=%d, i3=%d, fmaskpro均值=%f, himFT均值=%f\n', ...
                i, i3, mean(abs(fmaskpro(:))), mean(abs(himFT(:))));
        end
        if countimg == T % momentum method
            % 计算目标图像的动量（vobj = 动量系数×上一次动量 + 当前与历史状态差异）
            vobj = eta_obj.*vobj0 + (himFT - ObjT);
            himFT = ObjT + vobj; % 用动量更新目标图像
            
            vobj0 = vobj;                  
            ObjT = himFT;
            % 瞳孔函数的动量更新（同目标图像）
            vp = eta_p.*vp0 + (fmaskpro - PT);
            fmaskpro = PT + vp;
            
            vp0 = vp;
            PT = fmaskpro;
            
            countimg = 0;
        end
    end
    %him=ifft2(ifftshift(himFT));
    %figure;
    %set(gcf,'outerposition',get(0,'ScreenSize'))
    %subplot(121);imshow(abs(him(50:end-50,50:end-50)),[]);title('Amplitude');
    %subplot(122);imshow(angle(him(50:end-50,50:end-50)),[]);title('Phase');
end

him=ifft2(ifftshift(himFT));

end