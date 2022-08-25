function [y] = ADMM(D, s,alpha, lambda, rho, inverse_matrix,  MAX_ITER, TOL_DIST_X)
%% ADMM
    x0 = D'*s(:);
    y = x0;
    u = zeros(size(y));
    cnt = 1;
    distance = inf;

    while(cnt < MAX_ITER && distance > TOL_DIST_X) 
        cnt = cnt + 1;
        % solve the x subproblem
        x =  inverse_matrix*(x0 + rho*y- rho*u);
        
        % solve the y subproblem
        temp = soft_thr(x + u, lambda/rho);
        distance = norm(temp-y, 2);
        y = temp;
        
        % update u
        u = u + alpha*(x - y); 
        % update stopping criteria
    end
end