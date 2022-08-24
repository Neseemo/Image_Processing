function [angle, dist_from_origin] = polar_coord(P)
dist_from_origin = vecnorm(P);
angle = zeros(1,size(P,2));
angle(P(1,:)>0) = atan(P(2,P(1,:)>0)./P(1,P(1,:)>0));
angle(P(1,:)<=0) = atan(P(2,P(1,:)<=0)./P(1,P(1,:)<=0));
