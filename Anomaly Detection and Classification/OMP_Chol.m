function x = OMP_Chol(s, D, L, tau)
%lowtriang = dsp.LowerTriangularSolver;
%uptriang = dsp.UpperTriangularSolver;
[~, N] = size(D);
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
e = zeros(N,1);
while resNorm>tau && l-1<L
    %% SWEEP STEP: look for the column of D that matches at best noisySignal
    % compute the residual w.r.t. each column of D
    
    e(:) = resNorm^2 - (r'*D).^2; 
    % find the column of D that matches at best r, i.e. jStar = argmin(e(j))
    [~,jStar] = min(e);
    
    %% Initialization of R
    if l==1
        R(l,l) = 1;
    else
        %% Update the matrix u and R
        %opts.LT = true;
        u = linsolve(R(1:l-1,1:l-1),D_omega(:,1:l-1)'*D(:,jStar));
        %opts.LT = false;
        R(l,1:l) = [u', sqrt(1-u'*u)];
    end
    
    % update the support set
    omega(l) =  jStar;
    D_omega(:,l) = D(:,jStar);
    % update the coefficients by solving the triangular systems
    %opts.LT = true;
    temp = linsolve(R(1:l,1:l),y(omega(1:l)));
    %opts.LT = false;
    %opts.UT = true;
    x(omega(1:l)) = linsolve(R(1:l,1:l)',temp);
    %opts.UT = false;
    r = s - D*x;
    resNorm = norm(r, 2);
    l = l+1;
end
%% Remark 1
% It seems that the algorithm always runs for the maximum number of
% iteration possible, so a trick to further increase the code speed is to
% avoid computing the solutions of the 2 systems until the last iterations

%% REMARK 2
% I noticed that linsolve was faster when opts was not given as input. 
% This is probably because the systems are quite small and it's not
% necessary to tell the function that the system is triangular


