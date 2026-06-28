close all;clear;clc;

%% Prepare the experimental data
% Add necessary folders into the current working directory
addpath(genpath(pwd));
% Load data file
% data_name = 'MouseKidney_green';
% data_dir = 'FP_input_data_Sim.mat';

set(0, 'DefaultFigureVisible', 'on');

str_z = "z=_0.050_high_2nd";
ptr = "_up";
data_name = "FPM_RawData_z " + str_z + ptr;
batch_folder = "0619";
data_dir = "Raw_input\" + batch_folder + "\" + data_name + ".mat";

% Minimal external LED intensity correction configuration.
use_external_led_intensity_correction = true;
calib_batch_folder = "0605";
calib_file_name = "diffuser_calib_0605_毛玻璃01_5x5.mat";

load(data_dir); % refer to 'data_description.txt' for more details
imlow_HDR = im_high_HDR;
fprintf(' %s 的尺寸为 %d %d %d\n', data_name, size(imlow_HDR));

if use_external_led_intensity_correction
    calib_path = fullfile('Results', 'Calibration', calib_batch_folder, calib_file_name);
    calib_data = load(calib_path, 'calib');
    if ~isfield(calib_data, 'calib')
        error('Calibration file does not contain calib struct: %s', calib_path);
    end
    if ~isfield(calib_data.calib, 'led_intensity_fit_na_norm')
        error('Calibration file missing led_intensity_fit_na_norm: %s', calib_path);
    end

    led_intensity_rel = double(calib_data.calib.led_intensity_fit_na_norm(:));
    if numel(led_intensity_rel) ~= size(imlow_HDR, 3)
        error('LED calibration length (%d) does not match number of images (%d).', ...
            numel(led_intensity_rel), size(imlow_HDR, 3));
    end
    if any(led_intensity_rel <= 0) || any(~isfinite(led_intensity_rel))
        error('LED calibration contains non-positive or non-finite values.');
    end

    % imlow_HDR stores amplitude, so intensity calibration must be converted
    % to amplitude scaling: A_corr = A_raw / sqrt(I_rel).
    amp_scale = 1 ./ sqrt(led_intensity_rel);
    for k = 1:size(imlow_HDR, 3)
        imlow_HDR(:,:,k) = imlow_HDR(:,:,k) .* amp_scale(k);
    end
end

% Display raw images
is_show = 'center'; % 'center' shows the first low-res raw image; 'all' dynamically shows all low-res images
if strcmp(is_show,'center')
    figure();
    set(gcf,'outerposition',get(0,'ScreenSize'))
    imshow(imlow_HDR(:,:,1),[]);
    title(['raw image ' num2str(1)]);
% elseif strcmp(is_show,'all')
%     for k = 1:25%size(imlow_HDR,3)
%         figure(1);
%         set(gcf,'outerposition',get(0,'ScreenSize'))
%         imshow(imlow_HDR(:,:,k),[]);
%         title(['raw image ' num2str(k)]); pause(0.1);
%     end
end

%% Set up the experiment parameters
aberration = 0;
theta = 1.0;
xint = 0;
yint = 0;

xstart = 0; ystart = 0; % absolute coordinate of initial LED
arraysize = 5; % side length of lit LED array
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);
H      = 40.24; % distance between LEDs and sample, in mm
LEDp   = 4;     % distance between adjacent LEDs, in mm
nglass = 1.52;  % refraction index of glass substrate
t      = 1;     % glass thickness, in mm
[kx, ky, NAt] = k_vector(xlocation-xstart, ylocation-ystart, H, LEDp, nglass, t, theta, xint, yint, arraysize^2);

%% Reconstruct by FP algorithm
NA          = 0.15;     % objective NA
% spsize      = 1.38e-6;  % pixel size of low-res image on sample plane, in meter
upsmp_ratio = 4;        % upsampling ratio
psize       = spsize/upsmp_ratio; % pixel size of high-res image on sample plane, in meter

opts.loopnum    = 10;   % iteration number
opts.alpha      = 1;    % '1' for ePIE, other value for rPIE
opts.beta       = 1;    % '1' for ePIE, other value for rPIE
opts.gamma_obj  = 1;    % the step size for object updating
opts.gamma_p    = 1;    % the step size for pupil updating
opts.eta_obj    = 0.2;  % the step size for adding momentum to object updating
opts.eta_p      = 0.2;  % the step size for adding momentum to pupil updating
opts.T          = 1;    % do momentum every T images. '0' for no momentum during the recovery; integer, generally (0, arraysize^2].
opts.aberration = aberration; % pre-calibrated aberration, if available
opts.use_internal_intensity_correction = use_external_led_intensity_correction;

used_idx = 1:1:arraysize^2; % choose which raw image is used, for example, 1:2:arraysize^2 means do FPM recovery with No1 image, No3 image, No5 image......
imlow_used = imlow_HDR(:,:,used_idx);
kx_used = kx(used_idx);
ky_used = ky(used_idx);
[him, tt, fprobe, imlow_HDR1] = himrecover(imlow_used, kx_used, ky_used, NA, wlength, spsize, psize, z, opts);
% check result
%RMSE = show_recover(imlow_used, kx_used, ky_used, NA, wlength, spsize, psize, z, opts, him);
%disp(RMSE(1:2, 1:opts.loopnum + 2));
%
figure;
set(gcf,'outerposition',get(0,'ScreenSize'))
subplot(121);imshow(abs(him(50:end-50,50:end-50)),[]);title('Amplitude  High res');
subplot(122);imshow(angle(him(50:end-50,50:end-50)),[]);title('Phase  High res');

folder = fullfile('D:','1_徐前','物理','傅里叶叠层显微术','Code', 'Results', batch_folder);
if ~exist(folder, 'dir')
    mkdir(folder);
end
filename = "Recover_z " + str_z + ptr + ".png";
saveas(gcf, fullfile(folder, filename));

disp(['Wavelength: ',num2str(wlength.*1e+9),' nm, Loop: ',num2str(opts.loopnum)]);
disp(['Maximum illumination NA = ',num2str(max(NAt(used_idx)))]);

%% Save the results
% out_dir = "Results";
% mkdir(out_dir); addpath(out_dir);
% out_name = data_name + "_result" + ptr + ".mat";
% save(fullfile(folder, out_name), 'him','fprobe','tt','imlow_HDR1');
