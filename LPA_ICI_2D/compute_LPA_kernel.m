function g = compute_LPA_kernel(w, N)
% compute the LPA kernel for a given weights and polynomial degree
% input:
%   w: matrix containing the weights for the local LS problem (MxM patch)
%   N: degree of the polynomial approximation
% return:
%   g: the computed LPA kernel

% store the shape of the window function
r = size(w,1);
c = size(w, 2);

%% initialize the matrix T 

K = r*c;
% costruzione della matrice T
tx = linspace(0, 1, r)';
ty = linspace(0, 1, c)';
[tx, ty] = meshgrid(tx, ty);
tx = tx(:);
ty = ty(:);

T = zeros(K, (N+1)^2);
cnt = 1;
for lx=0:N
    for ly=0:(N - lx)
        if lx == 0 && ly == 0
            T(:, cnt) = ones(K,1);
        else
            T(:, cnt) = tx.^lx .* ty.^ly;
        end
        cnt = cnt + 1;
    end
end


%% use the code implemented in the previous lab!

% unrolled the weights
w = w(:);

% compute the kernel as in the 1D case
winv = zeros(size(w));
winv(w>0) = 1./w(w>0);
Winv = diag(winv);
W = diag(w);
[Q, ~] = qr(W*T, 0);
% define \tilde{Q}
Qtilde = Winv*Q;
% compute squared Qtilde
W2Qtilde =  W*W*Qtilde;
% select the central row of Qtilde
row = Qtilde(round(end/2),:);
g = zeros(K, 1);
for i = 1:min(N+1, K)
   g = g + W2Qtilde(:,i)*row(i); 
end

% before flipping, reshape the filter
g = reshape(g, r,c); 

% flipping
g = flip(g);

