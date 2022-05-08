% example of dct denoising
clear
close all

% load the image and rescale it in [0,1]
img = im2double(imread('cameraman.tif'));
% img = checkerboard(16, 8, 8);
% img = circshift(img, [1,1]);
% img = im2double(imread('BM3D_images/barbara.png'));
% img = im2double(imread('BM3D_images/Lena512.png'));

% noise level
sigma_noise = 20/255;

%patch size
p = 8;

% corrupt the image with white gaussian noise
% noisy_img = 

% compute the PSNR of the noisy input
% psnr_noisy = 

figure(1), imshow(img,[]), title('Original Image')
figure(2), imshow(noisy_img,[]), title(['Noisy Image, PSNR = ',num2str(psnr_noisy)])


%% generate the dct basis
D = zeros(p^2,p^2);

for ii=1:p^2
    patch = zeros(p,p);
    patch(ii) = 1;
    atom = idct2(patch);
    D(:,ii) = atom(:);
end

% display the 2D-DCT basis 
figure(), show_dictionary(D);
title(['atoms of 2D DCT ', num2str(p) ,' x ', num2str(p) , ' basis (unordered)'])

%% denoising

% initialize the estimated images
img_hat = zeros(size(noisy_img));

% initialize the weight matrix
weights = zeros(size(noisy_img));

% set the threshold for the Hart Thresholding
tau = 3 * sigma_noise; % Donoho says: sigma * sqrt(2*log(p^2))
 
% define the step (=p for non overlapping patches)
STEP = p; 

tic
% operates patch-wise
for ii = 1 : STEP : (size(noisy_img,1)-p+1)
    for jj=1 : STEP : (size(noisy_img,1)-p+1)
        % extrach the patch with the top left corner at pixel (ii, jj)
        % s  =
        
        % compute the representation w.r.t. the DCT dictionary
        % x =
        
        % perform Hard Thresholding (do not HT the dc component!)
        % x = 
        
        % perform the reconstruction
        % s_hat 
        
        % compute the weight for the reconstructed patch
        % w = 
        
        % put the reconstructed patch in the estimated image using the computed weight
        % UPDATE img_hat
        
        % store the weight of the current patch in the weight matrix
        % UPDATE weights
    end
end
toc

% normalize the estimated image with the computed weights
% img_hat =

% compute the psnr of the estimated image
% psnr_hat = 
str = sprintf('PSNR : %2.3f', psnr_hat);
disp(['\n', str, '\n'])
figure(3), imshow(img_hat,[]), title(['Estimated Image, ', str])

%% compare denoising by DCT agaist smoothing by convolution
% define a filter that by convolution returns the average over a 5x5 region
% h = 
% compute the convolution of the noisy image against the averaging filter
% img_hat_conv = 

% compute PSNR
% psnr_hat_conv =
% str = sprintf('PSNR : %2.3f', psnr_hat_conv);
% disp(['\n', str, '\n'])

% display the resulting image and check the loss of details
figure(4),
imshow(img_hat_conv,[]), title(['Estimated Image (conv), ', str])


%% noise estimation

% define the horizontal derivative filter
% h = 

% convolve the noisy_img and the filter

% compute sigma as the empirical std
% sigma_hat_emp =

% compute sigma using Median of Absolute Deviation
% sigma_hat =  
fprintf('sigma: %.3f, sigma_hat (empirical std): %.3f, sigma_hat (MAD): %.3f', sigma_noise, sigma_hat_emp, sigma_hat)



