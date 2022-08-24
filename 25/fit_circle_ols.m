function [theta] = fit_circle_ols(P)

%Given a Dataset P, finds through ols the best parameters "theta" that
%represent a circumference

%circumference eq: x^2 + y^2 + ax + by + c = 0
%theta = -[a,b,c]
A = P.^2; %x^2, y^2
A([3,4],:) = P; %x, y
A(5, :) = ones(1, size(P, 2)); %constant
theta = linsolve(A(3:end, :)*A(3:end, :)', A(3:end, :)*sum(A([1,2], :), 1)');

