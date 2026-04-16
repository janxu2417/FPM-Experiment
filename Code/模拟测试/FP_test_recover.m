function [] = FP_test_recover(imlow_HDR, kx, ky, NAt, spsize, wlength, NA_obj, z, opts, O)

% Display raw images
figure();
set(gcf,'outerposition',get(0,'ScreenSize'))
imshow(imlow_HDR(:,:,1),[]);
title(['raw image ' num2str(1)]);

%% Reconstruct by FP algorithm
%NA          = 0.1;      % objective NA
%spsize      = 1.845e-6; % pixel size of low-res image on sample plane, in meter
upsmp_ratio = 4;        % upsampling ratio
psize       = spsize/upsmp_ratio; % pixel size of high-res image on sample plane, in meter
arraysize = 10; % side length of lit LED array

used_idx = 1:1:arraysize^2; % choose which raw image is used, for example, 1:2:arraysize^2 means do FPM recovery with No1 image, No3 image, No5 image......
imlow_used = imlow_HDR(:,:,used_idx);
kx_used = kx(used_idx);
ky_used = ky(used_idx);
[him, tt, fprobe, imlow_HDR1] = himrecover(imlow_used, kx_used, ky_used, NA_obj, wlength, spsize, psize, z, opts);
% check result
%RMSE = show_recover(imlow_used, kx_used, ky_used, NA_obj, wlength, spsize, psize, z, opts, O);
%disp(RMSE(1:2, 1:opts.loopnum + 2));
%
figure;
set(gcf,'outerposition',get(0,'ScreenSize'))
subplot(121);imshow(abs(him(50:end-50,50:end-50)),[]);title('Amplitude  High res');
subplot(122);imshow(angle(him(50:end-50,50:end-50)),[]);title('Phase  High res');
disp(['Wavelength: ',num2str(wlength.*1e+9),' nm, Loop: ',num2str(opts.loopnum)]);
disp(['Maximum illumination NA = ',num2str(max(NAt(used_idx)))]);

ssim_abs = ssim(mat2gray(abs(him)),mat2gray(abs(O)));
ssim_angle = ssim(mat2gray(angle(him)),mat2gray(angle(O)));
disp(['振幅： ', num2str(ssim_abs)]);
disp(['相位： ', num2str(ssim_angle)]);

%% Save the results
%{
data_name = 'FP_rand';
out_dir = 'Results';
mkdir(out_dir); addpath(out_dir);
out_name = [data_name '_result.mat'];
save([out_dir,'\',out_name],'him','fprobe','tt','imlow_HDR1','ssim_abs','ssim_angle'); 
%}
end