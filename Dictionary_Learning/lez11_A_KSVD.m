clear 
close all
rng(1)
% load the image
img = im2double(imread('BM3D_images/barbara.png'));

imsz = size(img);

% patch size
p = 8;

% number of elements in the patch
M = p ^ 2;

% extract the random patches from the image
% use a small number (e.g., 1000) tu debug, then increase it
npatch = 1000;

% extract random patches from the image
S = zeros(M, npatch);
for i= 1:npatch
    ii = randsample(size(img, 1)-p+1, 1);
    jj = randsample(size(img,2)-p+1, 1);
    s = img(ii:ii+p-1,jj:jj+p-1);
    % bring all the patches to zero mean
    s = s - mean(s, "all");
    S(:,i) = reshape(s, M, 1);
end


% only few KSVD iterations are needed for a good dictionary
max_iter = 10;

% number of columns in the dictionary
N = 256;

% target sparsity
L = 4;

% initialize the dictionary randomly
D =  normrnd( 0 , 1 , 64, N);

% normalize the atoms of the dictionary
for i=1:N
    D(:,i) = D(:,i)/vecnorm(D(:,i));
end


% initialize the coefficient matrix
X = zeros(N, npatch);
ignored_counter = 0;

f = waitbar(0,'Please wait...');
for iter=1:max_iter
    if iter==1
        tic
    end
    % perform the sparse coding via OMP of all the columns of S
    for n=1:npatch
        X(:, n) = OMP_Chol(S(:,n), D, L, 0.001);  
        %X(:, n) = OMP(D, S(:,n), 0.001, L); 
    end
    % iterate over the columns of D
    for j=1:N

        % find which signals uses the j-th atom in the sparse coding
        omega = find(X(j,:) ~= 0);

           
        if isempty(omega)
            % if the atom is never used then ignore or substitute it
            ignored_counter = ignored_counter + 1;
            % compute the residual 
            s_hat = D * X;
            err = s_hat - S;
            res = vecnorm(err, 2);
            [~,omega] = max(abs(err));
        else
            % compute the residual matrix E, ignoring the j-th atom
            Dtemp = D;
            Dtemp(:,j) = zeros(p*p,1);
            E = S- Dtemp*X;

            % restrict E to the columns indicated by omega
            Eomega = E(:,omega);

            % compute the SVD of Eomega
            [U, Sigma, Vt] = svds(Eomega,1);

            % update the dictionary
            D(:,j) = U(:,1);

            % update the coefficient matrix
            X(j,omega) = Sigma(1,1) * Vt(1,:);

        end
    end
    if iter==1
        t=toc;
    end
    remaining_time = t*(max_iter-iter);
    waitbar(iter/max_iter,f,['Please wait ', num2str(remaining_time), ' seconds']);
end
delete(f)
toc
% show the learned dictionary
show_dictionary(D)



