function [bestmodel, bestinliers] = simpleMSAC(X, cardmss, epsi, do_debugPlot)
%SIMPLERANSAC - Robust fit with the MSAC algorithm
% adapted by Luca Magri from a code of Andrea Fusiello Computer Vision Toolkit

n = size(X,2);
alpha = 0.99; % Desired probability of success
f = 0.5 ; % Pessimistic estimate of inliers fraction

% set maximum number of iterations
MaxIterations  = log(1-alpha) / log(1 - (1 - f)^cardmss);
mincost =  +Inf;
A = X;
A(3, :) = ones(1, size(X, 2));

for  i = 1:MaxIterations
    
    % Generate s random indicies in the range 1..npts
    mss = randsample(n,cardmss); 
    % Fit model to this minimal sample set 
    [theta] = fit_line_dlt(X(:,mss));
    % Evaluate distances between points and model
    residuals = A'*theta;

    % Compute MSAC score
    idx = abs(residuals)<epsi;
    cost = sum(abs(residuals(idx)), "all") + sum(~ idx);
     
    % replace mincost, bestinliers and bestmodel if needed
    if cost < mincost
        mincost = cost;
        bestinliers = find(idx);
        bestmodel = theta;
    end
    
    % debug plots
    if(do_debugPlot)
        figure(99);
        clf;
        scatter(X(1,:),X(2,:));
        hold on;
        if cost <= mincost
            display_band(X, theta, epsi , 'g')
        else
            display_band(X, theta, epsi , 'r')
        end
        plot(X(1,mss), X(2, mss), 'kx', 'LineWidth', 3)
        xlim([0,1.2]);
        ylim([0.5,3.5]);
        title(['MSAC iterations: ', num2str(i), '/', num2str(MaxIterations), ' cost: ', num2str(cost)]);
        axis equal;
        pause(0.2)
    end
    
    
end
end

