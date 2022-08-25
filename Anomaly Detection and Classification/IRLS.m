function [x_current] = IRLS(D,DtD, s, lambda, MAX_ITER, TOL_DIST_X)
%% IRLS
    x0 = D'*s(:);
    x = x0;
    distanceX = inf;
    cnt = 1;
    delta = 1e-6;
    while(cnt < MAX_ITER && distanceX > TOL_DIST_X) 
        cnt = cnt + 1;
        w = diag(1./(abs(x) + delta));
        x_current = linsolve((DtD + lambda * w), x0); 
        %this update is computationally expensive
        % update stopping criteria
        distanceX = norm(x_current - x);
        x = x_current;
    end
end