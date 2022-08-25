# Anomaly Detection and Classification based on dictionary learning and sparse coding
references for the algorithms:
* ["Defect detection in SEM images of nanofibrous materials"](https://boracchi.faculty.polimi.it/docs/2017_Anomaly_Detection_SEM.pdf)
* ["Scale invariant anomaly detection with multiscale group sparse models"](https://www.researchgate.net/publication/307516129_Scale-invariant_anomaly_detection_with_multiscale_group-sparse_models)
* [“Robust face recognition via sparse representation”](https://ieeexplore.ieee.org/document/4483511)
The main files are "Lez20_A_AnomalyDetection" and "Lez20_B_Classification"

ell1 sparse coding:
1. ADMM: alternating direction method of multipliers sparse coding algorithm 
2. IRLS: iteratively reweighted least squares sparse coding algorithm
ell0 sparse coding:
1. MP: matching pursuit sparse coding
2. OMP: orthogonal matching pursuit
3. OMP_Chol: orthogonal matching pursuit taking advantage of choleski decomposition
dictionary learning algorithms:
1. ksvd: ksvd algorithm, can run with ell0 sparse coding algorithms
2. mod_dl: Method of optimal directions, can run together with ell1 sparse coding algorithms
utils:
1. show_dictionary: plots a given dictionary 
2. soft_thr: applies soft thresholding to a vector  (proximal mapping of l1 norm)

