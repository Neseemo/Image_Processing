function [g] = compute_LPA_kernel(w, N, ty)
    % generate the inverse of weights
    winv = zeros(size(w));
    winv(w>0) = 1./w(w>0);
    % set to zero weights that are Inf
    Winv = diag(winv);
    W = diag(w);
    
    M = size(w, 1);
    t = ty(1:M);
    T = zeros(size(t,1), N+1);
    for i=1:size(T,2)
        T(:,i) = t'.^(i-1);
    end
    
    [Q, ~] = qr(W*T, 0);

    % define \tilde{Q}
    Qtilde = Winv*Q;

    % compute squared Qtilde
    W2Qtilde =  W*W*Qtilde;

    % select the central row of Qtilde
    row = Qtilde(round(M/2),:);
    % calcolo il kernel
    g = zeros(M,1);
    for i = 1:min(N+1, M)
       g = g + W2Qtilde(:,i)*row(i); 
    end
    g = flip(g);
        
end