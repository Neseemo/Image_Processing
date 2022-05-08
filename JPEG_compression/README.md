# JPEG Compression
Jpeg compression algorithm consists in computing the representation of each patch of the image with respect to a particular orthonormal basis (I used DCT basis), "hard-threshold" the coefficients, keeping only the most important ones and finally reverse the transformation.
The DCT is separable so it allowes to speed up the computations.