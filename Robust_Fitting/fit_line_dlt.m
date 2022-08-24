function [theta, residuals, residual_error] = fit_line_dlt(P)
%
% P = [(x_1, y_1); ... ;(x_N, y_N)]
% points where the line y = mx+q should pass through
%
% Giacomo Boracchi
% 

% design matrix
P = P';
A = P;
A(:, 3) = ones(size(P, 1), 1);

% SVD 
[~, ~, V] = svd(A);

theta = V(:, end); 

residuals = A*theta;
residual_error = sum(residuals.^2);
