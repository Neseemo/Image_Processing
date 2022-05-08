# Image Denoising
In this notebook are presented two algorithms for image denoising:


*   Hard thresholding of DCT representation of the image patches 
*   smoothing by convolution (average of the pixels)

The first one works better but assumes the variance of the added noise to be known (wich is not obviously the case in reality). However, in the final part of the notebook a noise estimation technique is presented and, as shown, the estimated variance is very close to the actual one so the first algorith remains the best.
