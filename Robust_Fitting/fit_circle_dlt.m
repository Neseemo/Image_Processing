function [theta_circle, residuals, residual_error] = fit_circle_dlt(P)
%
% 1) Estimate the center of the circumference as the barycentre of the
% points
% 2) Convert in polar coordinates
% 3) Run dlt for lines in order to estimate: we expect to get a straight
% horizontal line since the radius should be constant. In this case the
% vertical distances from this line are also the residuals given by the
% "fit_line_dlt" method; since the straight line is assumed to horizontal,
% this are also the residuals in cartesian coordinates!
% 
%
% Remember eq of circumference: x^2 + y^2 + ax + by + c = 0
% theta = [a, b, c]. This functions returns -theta
estimated_center = mean(P,2);
[ang, d] = polar_coord(P - estimated_center);
P = [ang; d];
[theta, residuals, residual_error] = fit_line_dlt(P);
radius = theta(3)/theta(2); % the radius of the circle is the value of the line when x = 0
a = -2*estimated_center(1);
b = -2*estimated_center(2);
theta_circle = -[a, b, a^2/4+b^2/4-radius^2];
% return - theta 



