function x = OMP_Chol(s, D, L, tau)
%lowtriang = dsp.LowerTriangularSolver;
%uptriang = dsp.UpperTriangularSolver;
[M, N] = size(D);
x = zeros(N, 1);

% initialize the residual
r = s;

% initialize the residual norm
resNorm = norm(r,2);

% initialize the support set
omega = zeros(N,1);
D_omega = zeros(size(D));

% multiply the signal against D'
y = D'* s;
l = 1;
R = zeros(L,L);
while resNorm>tau && length(find(omega(1:l)~=0))<L
    %% SWEEP STEP: look for the column of D that matches at best noisySignal
    % compute the residual w.r.t. each column of D
    e = zeros(N,1);
    for j = 1 : N
        % For each direction dj compute the error
        % A large value is set for the direction that are already
        % in the active set in order to exclude them from the minimization
        % problem
        if any(omega(1:l) == j)
            e(j) = 99999; 
        else
            e(j) = resNorm^2 - (r' * D(:,j))^2;
        end
        
    end
    % find the column of D that matches at best r, i.e. jStar = argmin(e(j))
    [~,jStar] = min(e);
    
    %% Initialization of A and R
    if l==1
        % first iteration, initialize A
        %A = 1;
        R(l,l) = 1;
    else
        %% Update the matrix A
        opts.LT = true;
        %u = lowtriang(R(1:l-1,1:l-1),D_omega(:,1:l-1)'*D(:,jStar)); %v in the notes is (D_omega'*D(:,jStar)
        u = linsolve(R(1:l-1,1:l-1),D_omega(:,1:l-1)'*D(:,jStar), opts);
        opts.LT = false;
        R(l,1:l) = [u', sqrt(1-u'*u)];
        %release(lowtriang);
        %A = [A, u;u', 1];
    end
    
    % update the support set
    omega(l) =  jStar;
    D_omega(:,l) = D(:,jStar);
    % update the coefficients by solving the triangular systems
    opts.LT = true;
    %temp = lowtriang(R(1:l,1:l),y(omega(1:l)));
    temp = linsolve(R(1:l,1:l),y(omega(1:l)), opts);
    opts.LT = false;
    opts.UT = true;
    %x(omega(1:l)) = uptriang(R(1:l,1:l)',temp);
    x(omega(1:l)) = linsolve(R(1:l,1:l)',temp, opts);
    opts.UT = false;
    % remove the signal we have so far represented in coeff_MP (update the residual)
    r = s - D*x;
    resNorm = norm(r, 2);
    l = l+1;
    %rr=R(1:l-1,1:l-1)*R(1:l-1,1:l-1)'
    %dd=D_omega(:,1:l-1)'*D_omega(:,1:l-1)
    %release(lowtriang);
    %release(uptriang);
end
      


