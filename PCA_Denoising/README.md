# Image Denoising through PCA basis
The denoising task can be addressed also using PCA basis instead of DCT. 


*   The PRO is the more accurate representation we get since the principal components gather the most recurring patterns in the patches. 

*   The CON is the computational load required to compute the SVD decomposition (to deduce then the principal components). That's why when dealing with lots of images different methods are preferred.