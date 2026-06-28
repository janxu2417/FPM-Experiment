clear; close all; clc;

%% 彩色转成灰度图
out_dir = 'Data\实验数据';
out_name = 'expe6_result_2.mat';
%data_dir = 'Data\实验数据\expe5_result.mat';
%load(data_dir);
aberration = 0;
theta = 0;
xint = 0;
yint = 0;
z = 2.5000e-05;
wlength = 6.28e-07;
imlow_HDR = zeros(300, 300, 100);
end_j = 0;
for i = 15:2:214
    base_name = 'Data\实验数据\result_6_2\frame_';
    file_ext = '.png';
    %num_digits = 2;
    %读取彩色图片
    %num_str = sprintf(['%0', num2str(num_digits), 'd'], i);
    num_str = num2str(i);
    filename = [base_name, num_str, file_ext];
    rgbImage = imread(filename);

    %转换为灰度图
    grayImage = rgb2gray(rgbImage);

    end_j = end_j + 1;
    imlow_HDR(:,:,end_j) = double(grayImage(301:600,901:1200));
end
disp(['Wavelength: ',num2str(wlength.*1e+9),' nm, num: ',num2str(size(imlow_HDR,3))]);
imshow(imlow_HDR(:,:,1),[]);
save([out_dir,'\',out_name],'aberration', 'theta', 'xint', 'yint', 'z', 'wlength', 'imlow_HDR');
