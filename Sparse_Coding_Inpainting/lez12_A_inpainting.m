% example of dct denoising
clear
close all


% load the image and rescale it in [0,1]
img = im2double(imread('peppers256.png'));
figure(), imshow(img,[0, 1]), title('Original Image')

[nrow, ncol] = size(img);

% patch size
p = 8;
% input size
M = p^2;

% noise level
sigma_noise = 20 / 255;

% set the number of pixels to be removed at random
perc_of_removed_pixels = 0.25;

% corrupt the image with white gaussian noise
noisy_img = img + random('normal', 0, sigma_noise, size(img)); 

% arbitrarily remove pixels setting them to zero
% [modify  noisy_img]
r = randi([1, size(img,1)*size(img,2)],[floor(256*256*perc_of_removed_pixels), 1]);

% define a mask to keep track of which pixels have been removed
% msk = % a binary image of the same size of img having 0 in removed pixels
msk = ones(size(img));
msk(r) = 0;
figure(), imshow(noisy_img, [0, 1]), title(['Noisy Image before inpainting, PSNR = ',num2str(psnr(noisy_img,img))])
figure(), imshow(msk, [0, 1]), title(['Dead pixels'])
%% Use a redundant DCT dictionary learned from another images

data = load('dict_nat_img.mat');
D = data.D;
%add a column to keep track of the avg patch
new_col = ones(size(D,1), 1);
new_col = new_col / norm(new_col,2);
D = [D, new_col];
% add a constant atom to D, KSVD was trained over patches with zero mean


% update dictionary dimension
N = size(D, 1);

%figure(), show_dictionary(D), title('redundandt DCT dictionary')

%% SET stopping criteria of OMP
% orthogonal matching pursuit uses sparsity and errors as stopping criteria
L = 4;

% threshold on the error tau depends on the number of nonzero pixels

%% Inpainting based on OMP starts


% loop on image patches
step = 1;
img_hat = zeros(nrow, ncol);
weights = zeros(nrow, ncol);
f = waitbar(0,'Please wait...');
for ii=1:step:(nrow-p+1)
    if ii==1
        tic
    end
    for jj=1:step:(ncol-p+1)
        
        % patch extracted from the image
        s = noisy_img(ii:ii+p-1,jj:jj+p-1);
        
        
        % patch extracted from the mask
        m = msk(ii:ii+p-1,jj:jj+p-1);
        
        % design the projection operator over the current patch
        proj = diag(m(:));
        
        % tau should be proportional to the number of pixels remaining in the patch
        tau = 1.15 * sigma_noise *sqrt(sum(proj,"all")/M);
        
        % sparse coding w.r.t. PD the inpainted dictionary using L and tau as stopping criteria

        x = OMP(proj*D, s(:), tau, L);
        % reconstruction: synthesis w.r.t. D the dictionary yielding sparse representation
        s_hat = D*x;
        % keep track of all the estimates for aggregation 
        % UPDATE img_hat
        
        img_hat(ii:ii+p-1,jj:jj+p-1) = img_hat(ii:ii+p-1,jj:jj+p-1) + reshape(s_hat, [p,p]);
        % UPDATE weights

        weights(ii:ii+p-1,jj:jj+p-1) = weights(ii:ii+p-1,jj:jj+p-1) + ones(p,p);
    end
    if ii==1
        t=toc;
    end
    remaining_time = t*((nrow-p+1)-ii);
    waitbar(ii/(nrow-p+1),f,['Please wait ', num2str(remaining_time), ' seconds']);
end
delete(f)
% aggregation
img_hat = img_hat./weights;

%% show the result
corr_img = zeros(size(img));
corr_img(msk==1)=noisy_img(msk==1);
figure(), imshow(corr_img, [0, 1]), title(['corrupted img, PSNR = ',num2str(psnr(corr_img,img))])


figure(), imshow(img_hat), 
title(['Estimated Image, PSNR = ', num2str(psnr(img_hat,img))])

