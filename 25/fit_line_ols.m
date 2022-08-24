function [theta, residual_error] = fit_line_ols(P)
%
% P = [(x_1, y_1); ... ;(x_N, y_N)]
% points where the line y = mx+q should pass through
%
% Giacomo Boracchi
% 
% Least square solution of an overdetermined system 
% argmin ||A*x - y||^2  -> 
% \partial/x  (||Ax - y||^2) = 2 A'(Ax - y)
% solution by zeroing the dereivative: 2 A'(Ax - y) = 0
% x = (A'A)^(-1) A'y

% design matrix
P = P';
A = ones(size(P));
A(:, 1) = P(:, 1);

% straight line coefficient
y = P(:, 2);
theta = (A' * A) \ A' * y;


residuals = y - A*theta;
residual_error = sum(residuals.^2);
