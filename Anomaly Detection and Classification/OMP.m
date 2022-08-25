function [x_OMP_f, resNorm] = OMP(D, s, MINIMUM_RES_NORM, L)
    %% Orthogonal Matching Pursuit
    [M,N] = size(D);
    
    % initialization
    s_hat_OMP = zeros(M,1);
    r = s;
    l = 1;
    resNorm = norm(r,2);
    e = zeros(N,1);
    idx_omega_0 = true(N,1);
    idx_omega_1 = false(N,1);
    
    while resNorm>MINIMUM_RES_NORM && l<L + 1
        %% SWEEP STEP: look for the column of D that matches at best noisySignal
        % compute the residual w.r.t. each column of D (only for the
        % unexplored columns)
        
        %compute all the residuals in one-shot to speed up the computations
        e(idx_omega_0) = resNorm^2 - (r' * D(:,idx_omega_0)).^2;
        e(idx_omega_1) = inf;
        % find the column of D that matches at best r, i.e. jStar = argmin(e(j))
        [~,jStar] = min(e);
        
        % UPDATE the support set        
        idx_omega_0(jStar) = 0;
        idx_omega_1(jStar) = 1;
        
        % update the coefficients by solving the least square problem min ||D_omega x - s ||
        Dw1 = D(:,idx_omega_1);
        x_OMP = linsolve(Dw1'*Dw1,Dw1' * s);
        
        
        % remove the signal we have so far represented in x_OMP (update the residual)
        s_hat_OMP(:) = (Dw1) * x_OMP;
        r = s - s_hat_OMP;
        l = l + 1;
        resNorm = norm(r, 2);
    end
    x_OMP_f = zeros(size(D,2),1);
    x_OMP_f(idx_omega_1) = x_OMP;
    
end