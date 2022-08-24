function [bestmodel, bestinliers] = simpleLMEDS(X, cardmss, do_debugPlot)
%SIMPLELMS - Robust fit with the LMEDS algorithm
% adapted by Luca Magri from a code of Andrea Fusiello Computer Vision Toolkit

n = size(X,2);
alpha = 0.99; % Desired probability of success
f = 0.5; % Pessimistic estimate of inliers fraction


% set maximum number of iterations
MaxIterations  = log(1-alpha) / log(1 - (1 - f)^cardmss);
mincost = +Inf;
A = X;
A(3, :) = ones(1, size(X, 2));
for  i = 1:MaxIterations
    
    % Generate s random indicies in the range 1..npts
    mss = randsample(n,cardmss); 

    % Fit model to this minimal sample set.
    [theta] = fit_line_dlt(X(:,mss));
    % Evaluate distances between points and model
    residuals = A'*theta;
    % Compute LMS score
    cost = median(abs(residuals)); 
    
    % define inliner threshold (does make sense only when the model provides a good fit)
    % compute the standard deviation of distances (you can use MAD)
    scale = 1.4826*mad(abs(residuals),1);
    % instead of 3-sigma rule, we do 2.5-sigma rule.
    inliers = abs(residuals) < 2.5*scale;
    
    if cost < mincost
        mincost = cost;
        bestinliers = inliers;
        bestmodel = theta;
    end
    
    if(do_debugPlot)
        figure(101);
        clf;
        scatter(X(1,:),X(2,:));
        hold on;
        if cost <= mincost
            display_band(X, theta, 0.01 , 'g')
        else
            display_band(X, theta, 0.01 , 'r')
        end
        xlim([0,1.2]);
        ylim([0.5,3.5]);
        axis equal;
        title(['L-MEDS iterations: ', num2str(i), ' cost: ', num2str(cost)]);
        pause(0.2)
    end
end
end
