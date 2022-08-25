function D = ksvd(S, natoms, L, max_iter, algorithm)

    switch algorithm
        case "OMP"
            disp("Using Ksvd with OMP...")
        case "OMP_Chol"
            disp("Using Ksvd with OMP_Chol...")
        case "MP"
            disp("Using Ksvd with MP...")
        otherwise
            msg = "The parameter Algorithm must be one of the followings 'OMP', 'OMP_Chol', 'MP'";
            error(msg);
    end

    
    npatch = size(S, 2);
    % initialize the dictionary randomly
    D =  randn(size(S,1), natoms);

    % normalize the atoms of the dictionary
    D = D./vecnorm(D, 2);
    

    % initialize the coefficient matrix
    X = zeros(natoms, npatch);


    for iter=1:max_iter
        % perform the sparse coding via OMP of all the columns of S
        for n=1:npatch
            switch algorithm
                case "OMP"
                    X(:, n) = OMP(D, S(:,n), 0.001, L);
                case "OMP_Chol"
                    X(:, n) = OMP_Chol(S(:,n), D, L, 0.001);
                case "MP"
                    X(:, n) = MP(D, S(:,n), 0.001, L, 100);
                otherwise
                    msg = "The parameter Algorithm must be one of the followings 'OMP', 'OMP_Chol', 'MP'";
                    error(msg);
            end
        end
        % iterate over the columns of D
        for j=1:natoms
    
            % find which signals uses the j-th atom in the sparse coding
            omega = find(X(j,:) ~= 0);
    
               
            if isempty(omega)
                % if the atom is never used then ignore or substitute it
                % compute the residual 
                s_hat = D * X;
                err = s_hat - S;
                
                [~,new_col] = max(vecnorm(err, 2));
                % update the dictionary
                
                D(:,j) = S(:,new_col);
                omega = new_col;
            end
                % compute the residual matrix E, ignoring the j-th atom
                Dtemp = D;
                Dtemp(:,j) = 0;
                E = S- Dtemp*X;

                % restrict E to the columns indicated by omega
                Eomega = E(:,omega);

                % compute the SVD of Eomega
                [U, Sigma, Vt] = svd(Eomega, "econ");

                % update the dictionary
                D(:,j) = U(:,1);

                % update the coefficient matrix
                X(j,omega(1:size(Vt,2))) = Sigma(1,1) * Vt(1,:);
                X(j,omega(size(Vt,2):end)) = 0;
        end
    end
    
end