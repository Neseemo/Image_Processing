# LPA Intersection of Confidence Interval (for 1D signals)
This algorithm is designed to find the local polynomial approximation of a given signal but also optimizing the window dimension for each time step of the signal using "Intersection of Confidence Interval" Algorithm. This allows us to obtain a better estimate than using only LPA. The algorithm can also be improved using aggregation of the results given by running on the same signal LPA_ICI with asymmetric weights.

"lez22_A_LPA_ICI" and "lez22_B_LPA_ICI_aggregation" are the main live scripts
Functions:
1. compute_LPA_kernel: computes the LPA kernel given the weights, order of the polynomial and a time window
2. get_synthetic_signal: generate a random signal. There are different possible shapes
