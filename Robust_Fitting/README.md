# Robust Fitting
*   In this section is implemented RANSAC and some of its variants (Msac and Lmeds). This algorithms will be used to fit different models of the same "kind" (Multi-Model fitting)

*   The file demo_robustmmf_script.m contains the main script while demo_robustmmf.mlx contains the same file as a live editor.

*   You can run the code on different set of points and number of outliers in order to check the robustness.

Functions:
          1. addOutliersInBB:  adds the desired number of outliers
          2. display_band: plots the fitted linear model with the desired bandwidth
          3. display_clust: plots points highlighting the cluster they belong to (useful for multi model fitting to distinguish different models)
          4. fit_line_dlt: finds the parameters of the straight line that minimizes the distance between the points and the line (dlt)
          5. fit_line_ols: finds the parameters of the straight line that minimizes the sum of sq residuals (ols)
          6. fit_circle_dlt: finds the parameters of the circumference that minimizes the distance between the points and the circumference (dlt)
          7. fit_circle_ols: finds the parameters of the circumference through solving the over determined system through ols
          8. plot_circle: plots a circle given its parameters
          9. polar_coord: convert cartesian coord to polar coord
          10. res_line: compute the residual of a group of points with respect to a line
