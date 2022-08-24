function Y = addOutliersInBB(X,nOutliers, k)
%ADDOUTLIERSINBB generate nOutliers in the bounding box of the data
% the bounding box is expanded a little bit, the bigger the k the more it
% is extended
if(nargin< 3)
    k = 0.1;
end
xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));
wx = xmax-xmin;
wy = ymax -ymin;
dx = k*wx;
dy = k*wy;

Y = [X,[(xmax-xmin+2*dx ).*rand(1,nOutliers)+(xmin-dx); (ymax- ymin + 2*dy).*rand(1,nOutliers)+(ymin -dy)]];


