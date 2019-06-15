function[in, g, H] = gH_psd_inner(x, params)
    %gradient and hessian of log-det penalty for each clique
    n = length(x);
    X = smat(x);
    
    [L_X, p] = chol(X, 'lower');
    in_psd = ~p;
    in = in_psd;
    
    %Q = svecTransMat(size(X, 1));
    
    if in
        L_inv = inv(L_X);
        X_inv = L_inv'*L_inv;
        
        g = -svec(X_inv);
        
        %H = zeros(n, n);
        %for i = 1:n
            %y = sparse(n, 1);
            %y(i) = 1;
            
            %Y = smat(y);
            %H_curr = svec(X_inv * Y * X_inv);
            %H(i, :) = H_curr;
            
            
            %dominant computational expense. Take advantage of this.

                        
        %end
        Xkron = kron(X_inv, X_inv);
        H = params.Q * (Xkron * params.Q');
                
    else
        g = NaN;
        H = NaN;
        L = NaN;
    end
end