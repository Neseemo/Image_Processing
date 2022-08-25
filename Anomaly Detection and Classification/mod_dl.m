function D = mod_dl(S, natoms, regularization_parameter, max_iter, algorithm)


    npatch = size(S, 2);
    % initialize the dictionary randomly
    D =  randn(size(S,1), natoms);
    % normalize the atoms of the dictionary
    D = D./vecnorm(D, 2);
    % initialize the coefficient matrix
    X = zeros(natoms, npatch);
    %Algorithms parameters
    DtD = D'*D;
    MAX_ITER = 1e3;
    TOL_DIST_X = 1e-3;

    switch algorithm
        case "ADMM"
            rho = 1;
            inverse_matrix = 1/rho * eye(size(D,2)) - 1/rho*D'*((eye(size(D,1))*rho + D*D')\D);
            alpha = 1;
            disp("Using Mod with ADMM...")
        case "IRLS"
            disp("Using Mod with IRLS...")
        otherwise 
            msg = "The parameter Algorithm must be one of the followings 'ADMM', 'IRLS'";
            error(msg);
    end
    
    
    
    

    for iter = 1:max_iter
        %perform the sparse coding using any of the implemented algorithm for all the patches in S
        for n = 1:npatch
            s = S(:, n);
            switch algorithm
                case "ADMM"
                    x = ADMM(D, s,alpha, regularization_parameter, rho,inverse_matrix,  MAX_ITER, TOL_DIST_X);
                case "IRLS"
                    x = IRLS(D,DtD, s, regularization_parameter, MAX_ITER, TOL_DIST_X);
                otherwise 
                    msg = "The parameter Algorithm must be one of the followings 'ADMM', 'IRLS'";
                    error(msg);
            end

            X(:, n) = x;
        end
        %perform the MOD update
        D = S * X' /(X * X');
    
        %normalize the columns
        D = D./vecnorm(D,2);
    end
end