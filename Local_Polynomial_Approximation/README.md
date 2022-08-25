# Local Polynomial Approximation

This is a technique to fit a polynomial that approximates locally the behaviour of a function given some noisy observations in a certain window. 
The quality of the approximation depends both on the desired degree of the polynomial and on the window dimension (increasing the number of observations the variance of the prediction decreases).
Moreover it is also possible to weight differently the contribution of observations in the same window (Weighted LPA).

the files "lez21_A_LPA" and "lez21_B_WLPA" and the respective scripts are the main files to test the LPA and WLPA algorithms.
Functions:
1. compute_LPA_kernel: computes the LPA kernel given the weights, the order of the polynomial and the a time window.
2. get_syntethic_signal: generates randomly a signal. There are different options to generate "block", "parabolic", "cubic" signals



