function [bestmodel, bestinliers] = simpleRANSAC(X, cardmss, epsi, do_debugPlot)
%SIMPLERANSAC - Robust fit with the RANSAC algorithm
% adapted by Luca Magri from a code of Andrea Fusiello Computer Vision Toolkit
n = size(X,2);
alpha = 0.99; % Desired probability of success
f = 0.9; % Pessimistic estimate of outlier fraction

% set maximum number of iterations
MaxIterations  = log(1-alpha) / log(1 - (1 - f)^cardmss);
maxscore = -Inf; 
if cardmss == 2
    
    A = X;
    A(3, :) = ones(1, size(X, 2));
else
    A = X.^2; %x^2, y^2
    A([3,4],:) = X; %x, y
    A(5, :) = ones(1, size(X, 2)); %constant
end
    
i = 1;
while  i < MaxIterations
    % Generate s random indicies in the range 1..npts
    mss = randsample(n,cardmss);
    
    % Fit model to this minimal sample set.
    if cardmss == 2
        %fit lines
        [theta] = fit_line_dlt(X(:,mss));
        % Evaluate distances between points and model
        residuals = A'*theta;
        % identify inliers
        inliers = abs(residuals) < epsi;
        % assess contensus (the number of inliers)
        score = size(find(inliers), 1);
    end
    if cardmss == 3
        %fit circle
        theta = linsolve(A(3:end, mss)', sum(A([1,2], mss), 1)');
        center = [theta(1)/2;theta(2)/2];
        radius = sqrt(theta(1)^2/4 + theta(2)^2/4 + theta(3));
        dist_point_center = vecnorm(X - center, 2);
        dist_point_circle = abs(dist_point_center - radius);
        
        inliers = dist_point_circle < epsi;
        % assess contensus (the number of inliers)
        score = size(find(inliers), 2);
        
    end


    
    % replace maxscore, bestinliers and bestmodel if needed
    if score > maxscore
        maxscore = score;
        bestinliers = inliers;
        bestmodel = theta;
        f = (n - score)/n;
        MaxIterations = log(1-alpha) / log(1 - (1 - f)^cardmss);
    end
    
    % debug plots
    if(do_debugPlot)
        figure(99);
        clf;
        scatter(X(1,:),X(2,:), 8);
        hold on;
        if score >= maxscore
            if cardmss == 2
                display_band(X, theta, epsi , 'g')
            end
            if cardmss == 3
                plot_circle(center(1), center(2), radius, 'g');
            end
        else
            if cardmss == 2
                display_band(X, theta, epsi , 'r');
            end
            if cardmss == 3
                plot_circle(center(1), center(2), radius, 'r');
            end
        end
        plot(X(1,mss), X(2, mss), 'kx', 'LineWidth', 3)
        xlim([-5,5]);
        ylim([-5,5]);
        title(['Ransac iterations: ', num2str(i), ' consensus: ', num2str(score)]);
        axis equal;
        pause(0.2)
    end
    i = i+1;
    
    
end
disp(['Ransac Iterations: ', num2str(floor(MaxIterations))])

end

