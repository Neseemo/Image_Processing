function x = lez11_B_OMP_Cholesky(s, D, L, tau)

M, N = size(D);
x = zeros(N, 1);

% initialize the residual
r = s;

% initialize the residual norm
resNorm = norm(r,2);

% initialize the support set
omega = zeros(N,1);

% multiply the signal against D'
%y = 

while resNorm>tau && length(find(omega~=0))< L
    % sweep step
    %% SWEEP STEP: look for the column of D that matches at best noisySignal
    % compute the residual w.r.t. each column of D
    e = zeros(N,1);
    for j = 1 : N
        % For each direction dj compute the error
        % A large value is set for the direction that are already
        % in the active set in order to exclude them from the minimization
        % problem
        if omega(j) == 0
            e(j) = resNorm(l)^2 - (r' * D(:,j))^2;
        else
            e(j) = 9999;
        end
        
    end
    % find the column of D that matches at best r, i.e. jStar = argmin(e(j))
    [~,jStar] = min(e);
    
    
    if len(omega) == 0
        % first iteration, initialize A
        A = 1;
    else
        % update the matrix A
        R = chol(A);
        matrix = diag(omega);
        matrix = matrix(:,omega == 1);
        D_omega = D*matrix;
        u = R\(D_omega'*D(:,jStar));
        A = [A, u;u', 1];
    end
    
    % update the support set
    omega(jStar) = 1;
    
    % update the coefficients by solving the triangular systems
    x(omega) = (D_omega'*D_omega)\D_omega'*s;
    
    % remove the signal we have so far represented in coeff_MP (update the residual)
    r = s - D*x;
    resNorm = norm(r, 2);
    
end

        


