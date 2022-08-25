function [x_MP, resNorm] = MP(D, s, MINIMUM_RES_NORM, L, maxiter)
    %% Weak Matching Pursuit
    % initialization
    [~,N] = size(D);
    x_MP = zeros(N,1);
    r = s;
    l = 1;
    resNorm = norm(r, 2);
    e = zeros(N,1);

    % stoppint criteria: continue until the sparsity of the representation reaches L 
    %                    or as long as resNorm is above a minimum value
    %                    or as long as a maxium number of iterations have been reached
    
    while sum(x_MP ~= 0) < L && resNorm >MINIMUM_RES_NORM && l <maxiter
        % SWEEP STEP: look for the column of D that matches at best the
        % current residual
        
        %all the residual are computed in one-shot to speed up the computations
        e(:) = resNorm^2 - (r'*D).^2; 
        
        [~,jStar] = min(e);
        
        x_MP(jStar) = x_MP(jStar) + r'*D(:,jStar);
        
        % remove the signal we have so far represented in x_WMP (update the residual)
        s_hat_MP = D*x_MP;
        r = s - s_hat_MP;
        l = l + 1;
        
        % update the residual norm
        resNorm = norm(r, 2);
    end
end
